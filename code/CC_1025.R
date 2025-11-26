#!/usr/bin/env Rscript
# CC_pipeline_full_fixed.R
# Komplettes Pipeline-Script mit Feature-Progression (LASSO vs Ranger)
# und robuster Ranger-Prediction-Extraktion

suppressPackageStartupMessages({
  library(yaml)
  library(data.table)
  library(glmnet)
  library(ranger)
  library(pROC)
  library(foreach)
  library(doParallel)
  library(caret)
  library(lubridate)
  library(jsonlite)
  library(ggplot2)
  library(pheatmap)
})

options(bitmapType = "cairo") # PDF/Headless robust

# Helper: null-coalesce
`%||%` <- function(a, b) if (!is.null(a)) a else b

# -----------------------
# Helper functions
# -----------------------

read_config <- function(path) yaml::read_yaml(path)

make_output_dir <- function(base, run_name) {
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  out <- file.path(base, paste0(ts, "_", run_name))
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  dirs <- c("logs","models","figures","metrics","intermediates","configs")
  sapply(file.path(out, dirs), dir.create, recursive = TRUE, showWarnings = FALSE)
  out
}

find_first_file <- function(dir_path, pattern) {
  if (is.null(dir_path) || !dir.exists(dir_path)) stop("Directory not found: ", dir_path)
  files <- list.files(dir_path, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) stop("No file matching '", pattern, "' in ", dir_path)
  files[1]
}

detect_tsv_id_col <- function(dt, id_candidates) {
  for (c in id_candidates) if (c %in% colnames(dt)) return(c)
  colnames(dt)[1]
}

map_fam_pheno_to_binary <- function(pheno_vec) {
  pheno_vec[pheno_vec == -9] <- NA
  if (all(pheno_vec %in% c(1,2,NA))) {
    return(ifelse(pheno_vec == 2, 1, ifelse(pheno_vec == 1, 0, NA)))
  }
  if (all(pheno_vec %in% c(0,1,NA))) return(pheno_vec)
  ifelse(is.na(pheno_vec), NA, ifelse(pheno_vec > 1, 1, 0))
}

extract_pcs <- function(eig_dt, pcs_prefix = "PC", n_pcs_use = 10) {
  if (ncol(eig_dt) < 3) {
    warning("No PCs in eigenvec; returning ID only.")
    return(data.frame(ID = eig_dt[[2]]))
  }
  pc_count <- min(ncol(eig_dt) - 2, n_pcs_use)
  pc_names <- paste0(pcs_prefix, seq_len(pc_count))
  out <- data.frame(ID = as.character(eig_dt[[2]]))
  for (i in seq_len(pc_count)) out[[pc_names[i]]] <- as.numeric(eig_dt[[2 + i]])
  out
}

pc_outlier_mask <- function(pcs_df, sd_threshold = 6) {
  if (!all(c("PC1","PC2") %in% colnames(pcs_df))) return(rep(TRUE, nrow(pcs_df)))
  z1 <- as.numeric(scale(pcs_df$PC1))
  z2 <- as.numeric(scale(pcs_df$PC2))
  (abs(z1) <= sd_threshold) & (abs(z2) <= sd_threshold)
}

# Robust extractor for ranger prediction probabilities for class "1"
get_ranger_prob1 <- function(pred_obj, model = NULL) {
  preds <- pred_obj$predictions
  # if vector (regression or single-col), return numeric vector
  if (is.null(dim(preds))) return(as.numeric(preds))
  # if matrix with dimnames and column "1" exists, use it
  dn <- dimnames(preds)
  if (!is.null(dn) && !is.null(dn[[2]]) && "1" %in% dn[[2]]) return(as.numeric(preds[, "1"]))
  # try to infer from model levels (if model provided)
  if (!is.null(model) && !is.null(model$forest) && !is.null(model$forest$levels)) {
    lv <- as.character(model$forest$levels)
    if ("1" %in% lv) return(as.numeric(preds[, which(lv == "1")]))
    # fallback: if two levels, take second column as probability of positive class
    if (length(lv) == 2) return(as.numeric(preds[, 2]))
  }
  # final fallback: if >=2 columns, take second
  if (ncol(preds) >= 2) return(as.numeric(preds[, 2]))
  stop("Cannot extract class-1 probabilities from ranger predict output")
}

# -----------------------
# Main
# -----------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript CC_pipeline_full_fixed.R config.yaml")
cfg_file <- args[1]
cfg <- read_config(cfg_file)

out_dir_base <- cfg$output_dir %||% "outputs"
run_name <- cfg$run_name %||% "lasso_ranger_compare"
outdir <- make_output_dir(out_dir_base, run_name)
cat("Output dir:", outdir, "\n")

# Global containers
all_metrics <- list()
all_progression_results <- list()

# Parallel backend
ncores <- cfg$n_cores %||% 4
cl <- makeCluster(ncores)
registerDoParallel(cl)
cat("Registered parallel backend with", ncores, "workers\n")

# Iterate cohorts
for (cohort in cfg$cohorts) {
  name        <- cohort$name
  fam_dir     <- cohort$fam_eig_dir
  eig_dir     <- cohort$fam_eig_dir
  tsv_dir     <- cohort$tsv_dir
  tsv_pattern <- cohort$tsv_pattern %||% "\\.tsv$"

  cat(sprintf("\n--- Preparing cohort: %s ---\n", name))

  fam_file <- find_first_file(fam_dir, "\\.fam$")
  eig_file <- find_first_file(eig_dir, "\\.eigenvec$")
  tsv_file <- find_first_file(tsv_dir, tsv_pattern)

  fam <- fread(fam_file, header = FALSE, data.table = FALSE)
  eig <- fread(eig_file, header = FALSE, data.table = FALSE)
  tsv <- fread(tsv_file, header = TRUE, data.table = FALSE)

  fam[[2]] <- as.character(fam[[2]])
  eig[[2]] <- as.character(eig[[2]])

  tsv_id_col <- detect_tsv_id_col(tsv, cfg$id_candidates %||% c("ID","IID","Sample_ID"))
  if (tsv_id_col != "ID") {
    if ("ID" %in% colnames(tsv)) tsv$ID <- NULL
    setnames(tsv, tsv_id_col, "ID")
  }
  tsv$ID <- as.character(tsv$ID)

  # CASE_CONTROL
  case_col <- cfg$case_column %||% "CASE_CONTROL"
  if (!(case_col %in% colnames(tsv))) {
    if (ncol(fam) >= 6) {
      fam_pheno <- data.frame(ID = fam[[2]],
                              CASE_FROM_FAM = map_fam_pheno_to_binary(as.numeric(fam[[6]])),
                              stringsAsFactors = FALSE)
      tsv <- merge(tsv, fam_pheno, by = "ID", all.x = TRUE)
      tsv[[case_col]] <- tsv$CASE_FROM_FAM
      tsv$CASE_FROM_FAM <- NULL
      cat(sprintf("Info: '%s' added from .fam.\n", case_col))
    } else {
      stop("No phenotype in TSV and .fam lacks column 6.")
    }
  }

  # SEX
  if (!("SEX" %in% colnames(tsv))) {
    if (ncol(fam) >= 5) {
      sex_vec <- as.numeric(fam[[5]])
      tsv <- merge(tsv, data.frame(ID = fam[[2]], SEX = sex_vec), by = "ID", all.x = TRUE)
    } else {
      tsv$SEX <- NA_real_
      warning("SEX not found in .fam; set to NA.")
    }
  }
  tsv$SEX <- as.numeric(tsv$SEX)
  valid_sex_mask <- !is.na(tsv$SEX) & tsv$SEX %in% c(1, 2)
  removed_sex <- sum(!valid_sex_mask)
  if (removed_sex > 0) {
    cat(sprintf("SEX filtering: removed %d rows with SEX not in {1,2}.\n", removed_sex))
    tsv <- tsv[valid_sex_mask, , drop = FALSE]
  }

  # PCs
  pcs_prefix <- cfg$pc_prefix %||% "PC"
  n_pcs_use <- cfg$n_pcs_use %||% 10
  pcs_df <- extract_pcs(eig, pcs_prefix = pcs_prefix, n_pcs_use = n_pcs_use)
  tsv <- merge(tsv, pcs_df, by = "ID", all.x = TRUE)

  # Common IDs and ordering
  common_ids <- Reduce(intersect, list(fam[[2]], eig[[2]], tsv$ID))
  cat(sprintf("Common IDs: %d\n", length(common_ids)))
  if (length(common_ids) == 0) stop("No common IDs!")
  tsv <- tsv[tsv$ID %in% common_ids, , drop = FALSE]
  tsv <- tsv[order(tsv$ID), , drop = FALSE]

  # PC outlier removal
  pc1 <- paste0(pcs_prefix, 1); pc2 <- paste0(pcs_prefix, 2)
  if (all(c(pc1, pc2) %in% colnames(tsv))) {
    mask <- pc_outlier_mask(tsv[, c(pc1, pc2)], sd_threshold = cfg$outlier_sd %||% 6)
    removed <- sum(!mask)
    if (removed > 0) {
      cat(sprintf("PC outlier removal: %d samples removed (>%g SD).\n", removed, cfg$outlier_sd %||% 6))
      tsv <- tsv[mask, , drop = FALSE]
    }
  } else {
    cat("PC1/PC2 not present; skipping outlier removal.\n")
  }

  stopifnot(!any(duplicated(tsv$ID)))
  cat(sprintf("Final dataset for %s: %d samples x %d cols\n", name, nrow(tsv), ncol(tsv)))

  # Feature filter: keep only meta + Pt_1_ (or configured prefix)
  feature_prefix <- cfg$feature_prefix %||% "Pt_1_"
  feature_cols <- grep(paste0("^", feature_prefix), colnames(tsv), value = TRUE)
  meta_cols <- colnames(tsv)[!grepl("^Pt_", colnames(tsv))]
  tsv <- tsv[, c(meta_cols, feature_cols), drop = FALSE]
  cat("Kept", length(feature_cols), "features with prefix", feature_prefix, "and moved them to the end\n")

  # write prepared TSV for inspection
  out_file <- file.path(outdir, paste0(name, "_prepared.tsv"))
  fwrite(tsv, out_file, sep = "\t", quote = FALSE, na = "NA")
  cat("Written: ", out_file, "\n")

  # Prepare modelling inputs
  case_col <- cfg$case_column %||% "CASE_CONTROL"
  # ensure binary numeric 0/1 for general processing
  tsv[[case_col]] <- as.numeric(tsv[[case_col]])
  tsv <- tsv[!is.na(tsv[[case_col]]), , drop = FALSE]

  meta_keep <- c("ID", case_col, "SEX")
  meta_keep <- meta_keep[meta_keep %in% colnames(tsv)]
  feature_names <- setdiff(colnames(tsv), meta_keep)

  X_all <- as.data.frame(tsv[, feature_names, drop = FALSE])
  y_all <- as.numeric(tsv[[case_col]])

  cat("Data rows:", nrow(X_all), "columns:", ncol(X_all), "\n")

  # create outer folds
  n_outer <- cfg$n_outer %||% 3
  set.seed(cfg$seed %||% 42)
  outer_folds <- createFolds(y_all, k = n_outer, list = TRUE, returnTrain = FALSE)

  run_modes <- cfg$run_modes %||% list("combined")
  for (mode in run_modes) {
    cat("Running mode:", mode, "\n")
    subtag <- mode

    for (o in seq_along(outer_folds)) {
      test_idx_outer <- outer_folds[[o]]
      train_idx_outer <- setdiff(seq_len(nrow(X_all)), test_idx_outer)

      X_train_full <- X_all[train_idx_outer, , drop = FALSE]
      y_train_full <- y_all[train_idx_outer]
      X_test_full <- X_all[test_idx_outer, , drop = FALSE]
      y_test_full <- y_all[test_idx_outer]

      # inner folds
      inner_folds <- createFolds(y_train_full, k = cfg$n_inner %||% 3, list = TRUE, returnTrain = FALSE)

      for (fold_id in seq_along(inner_folds)) {
        val_idx <- inner_folds[[fold_id]]
        inner_train_idx <- setdiff(seq_len(nrow(X_train_full)), val_idx)

        X_train <- X_train_full[inner_train_idx, , drop = FALSE]
        y_train <- y_train_full[inner_train_idx]
        X_val <- X_train_full[val_idx, , drop = FALSE]
        y_val <- y_train_full[val_idx]

        X_test <- X_test_full
        y_test <- y_test_full

        # Baseline LASSO (all features)
        Xtr_mat <- as.matrix(X_train)
        Xte_mat <- as.matrix(X_test)
        train_center <- colMeans(Xtr_mat, na.rm = TRUE)
        train_scale <- apply(Xtr_mat, 2, sd, na.rm = TRUE)
        train_scale[train_scale == 0] <- 1
        Xtr_scaled <- scale(Xtr_mat, center = train_center, scale = train_scale)
        Xte_scaled <- scale(Xte_mat, center = train_center, scale = train_scale)

        cvfit_baseline <- cv.glmnet(x = Xtr_scaled, y = y_train, family = "binomial", alpha = 1, nfolds = 5, type.measure = "auc", parallel = FALSE)
        preds_lasso <- as.numeric(predict(cvfit_baseline, newx = Xte_scaled, s = "lambda.min", type = "response"))
        lasso_auc_baseline <- as.numeric(pROC::auc(y_test, preds_lasso))

        # Ranger baseline: ensure factor for classification
        y_train_factor <- factor(y_train, levels = c(0,1))
        rfit_baseline <- ranger::ranger(x = as.data.frame(Xtr_mat), y = y_train_factor, probability = TRUE, num.trees = cfg$num_trees %||% 500)
        pred_obj_baseline <- predict(rfit_baseline, data = as.data.frame(Xte_mat))
        preds_ranger <- get_ranger_prob1(pred_obj_baseline, rfit_baseline)
        ranger_auc_baseline <- as.numeric(pROC::auc(y_test, preds_ranger))

        lasso_metrics <- list(auc = lasso_auc_baseline, n_selected = NA)
        ranger_metrics <- list(auc = ranger_auc_baseline)

        cat(sprintf("Mode %s Outer %d Fold %d\n", mode, o, fold_id))
        cat("Correlation filter removed 0 features\n")
        cat(sprintf("LASSO AUC: %g selected: %s \n", lasso_metrics$auc, ifelse(is.null(lasso_metrics$n_selected), "NA", lasso_metrics$n_selected)))
        cat(sprintf("Ranger AUC: %g \n", ranger_metrics$auc))

        # Optional baseline ROC plot
        if (cfg$io$save_plots %||% TRUE) {
          roc_l <- pROC::roc(y_test, preds_lasso, quiet = TRUE)
          roc_r <- pROC::roc(y_test, preds_ranger, quiet = TRUE)
          p <- ggplot() +
            geom_line(aes(x = rev(roc_l$specificities), y = rev(roc_l$sensitivities)), color = "blue") +
            geom_line(aes(x = rev(roc_r$specificities), y = rev(roc_r$sensitivities)), color = "red") +
            labs(title = paste0("ROC ", subtag, " outer", o, " fold", fold_id),
                 x = "1 - Specificity", y = "Sensitivity") +
            theme_minimal() +
            annotate("text", x = 0.6, y = 0.2, label = paste0("LASSO AUC=", round(lasso_metrics$auc,3))) +
            annotate("text", x = 0.6, y = 0.1, label = paste0("Ranger AUC=", round(ranger_metrics$auc,3)))
          ggsave(filename = file.path(outdir, "figures", paste0(subtag, "_outer", o, "_fold", fold_id, "_roc_baseline.pdf")),
                 plot = p, width = 6, height = 4)
        }

        # --- Feature-Progression Block (LASSO vs Ranger) ---
        fea_list <- cfg$fea_list %||% c(1,3,5,7,9,10,12,15,20,25,30)
        fea_list <- sort(unique(as.integer(fea_list[fea_list > 0])))

        results_list <- data.table(method = character(), n_features = integer(), auc = numeric(), outer = integer(), fold = integer())

        # Prepare scaled matrices for LASSO ranking
        Xtr_mat <- as.matrix(X_train)
        Xte_mat <- as.matrix(X_test)
        train_center <- colMeans(Xtr_mat, na.rm = TRUE)
        train_scale <- apply(Xtr_mat, 2, sd, na.rm = TRUE)
        train_scale[train_scale == 0] <- 1
        Xtr_scaled <- scale(Xtr_mat, center = train_center, scale = train_scale)
        Xte_scaled <- scale(Xte_mat, center = train_center, scale = train_scale)

        # LASSO ranking
        cvfit <- cv.glmnet(x = Xtr_scaled, y = y_train, family = "binomial", alpha = 1, nfolds = 5, type.measure = "auc", parallel = FALSE)
        lasso_coefs <- as.numeric(coef(cvfit, s = "lambda.min")[-1])
        names(lasso_coefs) <- colnames(Xtr_mat)
        lasso_rank <- names(sort(abs(lasso_coefs), decreasing = TRUE))

        # Ranger ranking (train with factor target to ensure classification importance is computed)
        y_train_factor_full <- factor(y_train, levels = c(0,1))
        imp_method <- cfg$ranger$importance %||% "impurity"
        rfit_full <- ranger::ranger(x = as.data.frame(Xtr_mat), y = y_train_factor_full, probability = TRUE, importance = imp_method, num.trees = cfg$num_trees %||% 500)
        ranger_imp <- rfit_full$variable.importance
        ranger_rank <- names(sort(ranger_imp, decreasing = TRUE))

        for (k in fea_list) {
          k_use <- min(k, ncol(Xtr_mat))
          if (k_use <= 0) next

          # LASSO on top-k
          sel_lasso <- intersect(lasso_rank, colnames(Xtr_mat))[seq_len(k_use)]
          if (length(sel_lasso) > 0) {
            Xtr_lasso_k <- Xtr_scaled[, sel_lasso, drop = FALSE]
            Xte_lasso_k <- Xte_scaled[, sel_lasso, drop = FALSE]
            fit_k <- cv.glmnet(x = Xtr_lasso_k, y = y_train, family = "binomial", alpha = 1, nfolds = 5, type.measure = "auc", parallel = FALSE)
            preds_lasso_k <- as.numeric(predict(fit_k, newx = Xte_lasso_k, s = "lambda.min", type = "response"))
            auc_lasso_k <- as.numeric(pROC::auc(y_test, preds_lasso_k))
            results_list <- rbind(results_list, data.table(method = "LASSO", n_features = k_use, auc = auc_lasso_k, outer = o, fold = fold_id))
          }

          # Ranger on top-k
          sel_ranger <- intersect(ranger_rank, colnames(Xtr_mat))[seq_len(k_use)]
          if (length(sel_ranger) > 0) {
            Xtr_r_k <- as.data.frame(Xtr_mat[, sel_ranger, drop = FALSE])
            Xte_r_k <- as.data.frame(Xte_mat[, sel_ranger, drop = FALSE])
            # train with factor
            rfit_k <- ranger::ranger(x = Xtr_r_k, y = factor(y_train, levels = c(0,1)), probability = TRUE, importance = "none", num.trees = cfg$num_trees %||% 500)
            pred_obj_k <- predict(rfit_k, data = Xte_r_k)
            preds_ranger_k <- get_ranger_prob1(pred_obj_k, rfit_k)
            auc_ranger_k <- as.numeric(pROC::auc(y_test, preds_ranger_k))
            results_list <- rbind(results_list, data.table(method = "Ranger", n_features = k_use, auc = auc_ranger_k, outer = o, fold = fold_id))
          }
        } # end k loop

        # store fold progression
        keyname <- paste0(name, "_", subtag, "_outer", o, "_fold", fold_id)
        all_progression_results[[keyname]] <- results_list

        # plot AUC vs features per fold
        if ((cfg$io$save_plots %||% TRUE) && nrow(results_list) > 0) {
          p_prog <- ggplot(results_list, aes(x = n_features, y = auc, color = method)) +
            geom_line(aes(group = method)) +
            geom_point() +
            scale_x_continuous(breaks = fea_list) +
            theme_minimal() +
            labs(title = paste0("AUC vs #features: ", subtag, " outer", o, " fold", fold_id),
                 x = "Number of top features", y = "AUC")
          ggsave(filename = file.path(outdir, "figures", paste0(subtag, "_outer", o, "_fold", fold_id, "_auc_vs_features.pdf")),
                 plot = p_prog, width = 6, height = 4)
        }

        # summarize best per method and append to all_metrics
        best_per_method <- results_list[, .SD[which.max(auc)], by = method]
        l_entry <- list(
          lasso = list(method = "LASSO", subtag = subtag, outer = o, fold = fold_id,
                       auc = best_per_method[method == "LASSO"]$auc %||% lasso_metrics$auc,
                       auc_ci = c(NA, NA), time_sec = NA, n_selected = best_per_method[method == "LASSO"]$n_features %||% NA),
          ranger = list(method = "Ranger", subtag = subtag, outer = o, fold = fold_id,
                        auc = best_per_method[method == "Ranger"]$auc %||% ranger_metrics$auc,
                        auc_ci = c(NA, NA), time_sec = NA)
        )
        all_metrics <- c(all_metrics, list(l_entry))

      } # inner folds
    } # outer folds
  } # modes
} # cohorts

stopCluster(cl)

# Aggregate metrics
metrics_dt <- rbindlist(lapply(all_metrics, function(x) {
  rbindlist(list(
    data.table(method = x$lasso$method, subtag = x$lasso$subtag, outer = x$lasso$outer, fold = x$lasso$fold,
               auc = x$lasso$auc, auc_ci_lower = x$lasso$auc_ci[1], auc_ci_upper = x$lasso$auc_ci[2], time = x$lasso$time_sec,
               n_selected = x$lasso$n_selected),
    data.table(method = x$ranger$method, subtag = x$ranger$subtag, outer = x$ranger$outer, fold = x$ranger$fold,
               auc = x$ranger$auc, auc_ci_lower = x$ranger$auc_ci[1], auc_ci_upper = x$ranger$auc_ci[2], time = x$ranger$time_sec,
               n_selected = NA)
  ))
}))

fwrite(metrics_dt, file.path(outdir, "metrics", "aggregate_metrics.csv"))

# Summary plot
p1 <- ggplot(metrics_dt, aes(x = method, y = auc, fill = method)) +
  geom_boxplot() + theme_minimal() + ggtitle("AUC distribution per method")
ggsave(filename = file.path(outdir, "figures", "auc_distribution_methods.pdf"), plot = p1, width = 6, height = 4)

# Save progressions and metrics
saveRDS(all_progression_results, file.path(outdir, "metrics", "all_progression_results.rds"))
saveRDS(metrics_dt, file.path(outdir, "metrics", "metrics_dt.rds"))

cat("Pipeline finished. Results in:", outdir, "\n")
