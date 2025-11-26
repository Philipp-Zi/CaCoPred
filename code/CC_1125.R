#!/usr/bin/env Rscript
# CC_pipeline_full_unabridged.R
# Vollst�ndige wissenschaftliche Pipeline
# - Korrelation-basierte Feature-Selektion
# - LASSO vs Ranger (Ranger nur als finaler Pr�diktor, importance optional nur f�r Ranking)
# - Cross-cohort, sex stratify, robust handling of k=1, periodic flush to disk
# - Aggregated AUC vs #features plots (no per-fold plots by default)

suppressPackageStartupMessages({
  library(yaml)
  library(data.table)
  library(glmnet)
  library(ranger)
  library(pROC)
  library(foreach)
  library(doParallel)
  library(caret)
  library(ggplot2)
})

options(bitmapType = "cairo")

`%||%` <- function(a,b) if(!is.null(a)) a else b

# ----------------------- Helpers -----------------------

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
  if (all(pheno_vec %in% c(1,2,NA))) return(ifelse(pheno_vec == 2, 1, ifelse(pheno_vec == 1, 0, NA)))
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
  z1 <- as.numeric(scale(pcs_df$PC1)); z2 <- as.numeric(scale(pcs_df$PC2))
  (abs(z1) <= sd_threshold) & (abs(z2) <= sd_threshold)
}

# robust extractor for ranger pred proba for class "1"
get_ranger_prob1 <- function(pred_obj, model = NULL) {
  preds <- pred_obj$predictions
  if (is.null(dim(preds))) return(as.numeric(preds))
  dn <- dimnames(preds)
  if (!is.null(dn) && !is.null(dn[[2]]) && "1" %in% dn[[2]]) return(as.numeric(preds[, "1"]))
  if (!is.null(model) && !is.null(model$forest) && !is.null(model$forest$levels)) {
    lv <- as.character(model$forest$levels)
    if ("1" %in% lv) return(as.numeric(preds[, which(lv == "1")]))
    if (length(lv) == 2) return(as.numeric(preds[, 2]))
  }
  if (ncol(preds) >= 2) return(as.numeric(preds[, 2]))
  stop("Cannot extract class-1 probabilities from ranger predict output")
}

# safe LASSO fit/predict (handles 0/1/>=2 features)
fit_predict_lasso_safe <- function(Xtr, ytr, Xte, yte, nfolds = 5, alpha = 1) {
  # Remove constant columns
  nzv <- apply(Xtr, 2, function(z) sd(z, na.rm = TRUE) > 0)
  if (sum(nzv) == 0) return(list(preds = rep(NA_real_, length(yte)), auc = NA_real_))
  Xtr2 <- Xtr[, nzv, drop = FALSE]; Xte2 <- Xte[, nzv, drop = FALSE]
  if (ncol(Xtr2) == 1) {
    df <- data.frame(y = ytr, x = as.numeric(Xtr2[,1]))
    glm_fit <- stats::glm(y ~ x, data = df, family = binomial(link = "logit"))
    newdf <- data.frame(x = as.numeric(Xte2[,1]))
    preds <- stats::predict(glm_fit, newdata = newdf, type = "response")
    auc <- tryCatch(as.numeric(pROC::auc(yte, preds)), error = function(e) NA_real_)
    return(list(preds = preds, auc = auc))
  }
  fit <- cv.glmnet(x = Xtr2, y = ytr, family = "binomial", alpha = alpha, nfolds = nfolds, type.measure = "auc", parallel = FALSE)
  preds <- as.numeric(predict(fit, newx = Xte2, s = "lambda.min", type = "response"))
  auc <- tryCatch(as.numeric(pROC::auc(yte, preds)), error = function(e) NA_real_)
  return(list(preds = preds, auc = auc, model = fit))
}

# correlation-based feature filter
correlation_feature_filter <- function(X, threshold = 0.95, method = "pearson") {
  # X: numeric matrix or data.frame of candidate features
  # Remove columns with too many NAs or zero variance before correlation computation is caller's job
  if (ncol(X) <= 1) return(colnames(X))
  # compute correlation matrix (pairwise complete)
  cm <- abs(stats::cor(X, use = "pairwise.complete.obs", method = method))
  diag(cm) <- 0
  keep <- rep(TRUE, ncol(X))
  names(keep) <- colnames(X)
  # greedy remove: for each pair with corr > threshold remove the variable with larger mean abs corr
  while (TRUE) {
    high <- which(cm > threshold, arr.ind = TRUE)
    if (nrow(high) == 0) break
    # compute mean absolute corr per variable
    mac <- rowMeans(cm, na.rm = TRUE)
    # find pair to remove: choose variable with larger mac among first pair
    i <- high[1,1]; j <- high[1,2]
    if (mac[i] >= mac[j]) {
      rm_idx <- i
    } else rm_idx <- j
    # set corresponding column/row to zero and mark keep false
    keep[rm_idx] <- FALSE
    cm[rm_idx, ] <- 0; cm[, rm_idx] <- 0
  }
  names(keep)[keep]
}

# periodic flush helper: append dt to csv and clear dt
flush_progression_to_disk <- function(dt, outdir, part_idx) {
  f <- file.path(outdir, "metrics", paste0("progression_partial_", part_idx, ".csv"))
  fwrite(dt, f, append = FALSE)
}

# ----------------------- Main -----------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript CC_pipeline_full_unabridged.R config.yaml")
cfg_file <- args[1]
cfg <- read_config(cfg_file)

out_dir_base <- cfg$output_dir %||% cfg$io$output_base %||% "outputs"
run_name <- cfg$run_name %||% "lasso_ranger_compare"
outdir <- make_output_dir(out_dir_base, run_name)
cat("Output dir:", outdir, "\n")

# parallel backend
ncores <- cfg$n_cores %||% cfg$parallel$n_workers %||% 4
cl <- makeCluster(ncores); registerDoParallel(cl)
cat("Registered parallel backend with", ncores, "workers\n")

# containers
# all_progression_dt holds per-fold and cross eval progression rows
all_progression_dt <- data.table(cohort = character(), mode = character(), outer = integer(), fold = integer(), method = character(), n_features = integer(), auc = numeric())
flush_part_idx <- 0
flush_threshold <- cfg$io$flush_rows_threshold %||% 5000

# Load and prepare each cohort (merge fam/eigenvec/tsv, keep features)
cohort_data_list <- list()
for (cohort in cfg$cohorts) {
  name <- cohort$name
  cat("\n--- Preparing cohort:", name, "---\n")
  fam_file <- find_first_file(cohort$fam_eig_dir, "\\.fam$")
  eig_file <- find_first_file(cohort$fam_eig_dir, "\\.eigenvec$")
  tsv_file <- find_first_file(cohort$tsv_dir, cohort$tsv_pattern %||% "\\.ts$")
  fam <- fread(fam_file, header = FALSE, data.table = FALSE); eig <- fread(eig_file, header = FALSE, data.table = FALSE)
  tsv <- fread(tsv_file, header = TRUE, data.table = FALSE)
  fam[[2]] <- as.character(fam[[2]]); eig[[2]] <- as.character(eig[[2]])

  tsv_id_col <- detect_tsv_id_col(tsv, cfg$id_candidates %||% c("ID","IID","Sample_ID"))
  if (tsv_id_col != "ID") {
    if ("ID" %in% colnames(tsv)) tsv$ID <- NULL
    setnames(tsv, tsv_id_col, "ID")
  }
  tsv$ID <- as.character(tsv$ID)

  case_col <- cfg$case_column %||% "CASE_CONTROL"
  if (!(case_col %in% colnames(tsv))) {
    if (ncol(fam) >= 6) {
      fam_pheno <- data.frame(ID = fam[[2]], CASE_FROM_FAM = map_fam_pheno_to_binary(as.numeric(fam[[6]])), stringsAsFactors = FALSE)
      tsv <- merge(tsv, fam_pheno, by = "ID", all.x = TRUE)
      tsv[[case_col]] <- tsv$CASE_FROM_FAM; tsv$CASE_FROM_FAM <- NULL
      cat("Info: CASE_CONTROL added from .fam\n")
    } else stop("No phenotype in TSV and .fam lacks column 6.")
  }

  # SEX
  if (!("SEX" %in% colnames(tsv))) {
    if (ncol(fam) >= 5) {
      sex_vec <- as.numeric(fam[[5]])
      tsv <- merge(tsv, data.frame(ID = fam[[2]], SEX = sex_vec), by = "ID", all.x = TRUE)
    } else tsv$SEX <- NA_real_
  }
  tsv$SEX <- as.numeric(tsv$SEX)
  # PCs
  pcs_prefix <- cfg$pc_prefix %||% cfg$preproc$pc_prefix %||% "PC"
  n_pcs_use <- cfg$n_pcs_use %||% 10
  eig_pcs <- extract_pcs(eig, pcs_prefix = pcs_prefix, n_pcs_use = n_pcs_use)
  tsv <- merge(tsv, eig_pcs, by = "ID", all.x = TRUE)

  # common IDs
  common_ids <- Reduce(intersect, list(fam[[2]], eig[[2]], tsv$ID))
  cat("Common IDs:", length(common_ids), "\n")
  tsv <- tsv[tsv$ID %in% common_ids, , drop = FALSE]
  tsv <- tsv[order(tsv$ID), , drop = FALSE]

  # optional PC outlier removal (uses PC1/PC2)
  pc1 <- paste0(pcs_prefix, 1); pc2 <- paste0(pcs_prefix, 2)
  if (all(c(pc1, pc2) %in% colnames(tsv))) {
    mask <- pc_outlier_mask(tsv[, c(pc1, pc2)], sd_threshold = cfg$outlier_sd %||% 6)
    removed <- sum(!mask)
    if (removed > 0) {
      cat(sprintf("PC outlier removal: %d samples removed\n", removed))
      tsv <- tsv[mask, , drop = FALSE]
    }
  }

  # feature selection by prefix
  feature_prefix <- cfg$feature_prefix %||% cfg$preproc$feature_prefix %||% "Pt_1_"
  feature_cols <- grep(paste0("^", feature_prefix), colnames(tsv), value = TRUE)
  meta_cols <- colnames(tsv)[!grepl(paste0("^", feature_prefix), colnames(tsv))]
  tsv_prepared <- tsv[, c(meta_cols, feature_cols), drop = FALSE]
  fwrite(tsv_prepared, file.path(outdir, paste0(name, "_prepared.tsv")), sep = "\t", quote = FALSE, na = "NA")
  cat("Written prepared tsv for cohort:", name, "\n")
  cohort_data_list[[name]] <- list(tsv = tsv_prepared, feature_cols = feature_cols)
}

# ----------------------- Core pipeline -----------------------

# iterate each cohort (within-cohort cv + feature progression)
for (cohort_name in names(cohort_data_list)) {
  data_obj <- cohort_data_list[[cohort_name]]$tsv
  feature_cols <- cohort_data_list[[cohort_name]]$feature_cols
  case_col <- cfg$case_column %||% "CASE_CONTROL"
  meta_keep <- c("ID", case_col, "SEX")
  meta_keep <- meta_keep[meta_keep %in% colnames(data_obj)]
  feature_names <- setdiff(colnames(data_obj), meta_keep)
  X_all <- as.data.frame(data_obj[, feature_names, drop = FALSE])
  y_all <- as.numeric(data_obj[[case_col]])
  cat(sprintf("\n--- Cohort %s: %d samples, %d features ---\n", cohort_name, nrow(X_all), ncol(X_all)))

  # outer folds
  n_outer <- cfg$n_outer %||% cfg$cv$n_outer %||% 3
  set.seed(cfg$seed %||% 42)
  outer_folds <- createFolds(y_all, k = n_outer, list = TRUE, returnTrain = FALSE)

  # run modes (combined or sex strats)
  run_modes <- if (isTRUE(cfg$sex_stratify)) c("combined","male","female") else c("combined")

  for (mode in run_modes) {
    for (o in seq_along(outer_folds)) {
      test_idx_outer <- outer_folds[[o]]; train_idx_outer <- setdiff(seq_len(nrow(X_all)), test_idx_outer)
      X_train_full <- X_all[train_idx_outer, , drop = FALSE]; y_train_full <- y_all[train_idx_outer]
      X_test_full <- X_all[test_idx_outer, , drop = FALSE]; y_test_full <- y_all[test_idx_outer]

      # apply sex stratification if requested
      if (mode == "male") {
        sel_train <- which(data_obj$SEX[train_idx_outer] == 1); sel_test <- which(data_obj$SEX[test_idx_outer] == 1)
        if (length(sel_train) < (cfg$development$min_strata_size %||% 10) || length(sel_test) < (cfg$development$min_strata_size %||% 10)) { cat("Skipping male due to small N\n"); next }
        X_train_full <- X_train_full[sel_train, , drop = FALSE]; y_train_full <- y_train_full[sel_train]
        X_test_full <- X_test_full[sel_test, , drop = FALSE]; y_test_full <- y_test_full[sel_test]
      } else if (mode == "female") {
        sel_train <- which(data_obj$SEX[train_idx_outer] == 2); sel_test <- which(data_obj$SEX[test_idx_outer] == 2)
        if (length(sel_train) < (cfg$development$min_strata_size %||% 10) || length(sel_test) < (cfg$development$min_strata_size %||% 10)) { cat("Skipping female due to small N\n"); next }
        X_train_full <- X_train_full[sel_train, , drop = FALSE]; y_train_full <- y_train_full[sel_train]
        X_test_full <- X_test_full[sel_test, , drop = FALSE]; y_test_full <- y_test_full[sel_test]
      }

      # inner folds for nested CV (kept as before)
      inner_k <- cfg$n_inner %||% cfg$cv$n_inner %||% 3
      inner_folds <- createFolds(y_train_full, k = inner_k, list = TRUE, returnTrain = FALSE)

      for (fold_id in seq_along(inner_folds)) {
        val_idx <- inner_folds[[fold_id]]; inner_train_idx <- setdiff(seq_len(nrow(X_train_full)), val_idx)
        X_train <- X_train_full[inner_train_idx, , drop = FALSE]; y_train <- y_train_full[inner_train_idx]
        X_val <- X_train_full[val_idx, , drop = FALSE]; y_val <- y_train_full[val_idx]
        X_test <- X_test_full; y_test <- y_test_full

        # Standardize using train
        Xtr_mat <- as.matrix(X_train); Xte_mat <- as.matrix(X_test)
        train_center <- colMeans(Xtr_mat, na.rm = TRUE)
        train_scale <- apply(Xtr_mat, 2, sd, na.rm = TRUE); train_scale[train_scale == 0] <- 1
        Xtr_scaled <- scale(Xtr_mat, center = train_center, scale = train_scale)
        Xte_scaled <- scale(Xte_mat, center = train_center, scale = train_scale)

        # Pre-feature-filter: correlation-based filter applied on inner-train (scientifically correct)
        corr_threshold <- cfg$preproc$correlation_filter_threshold %||% cfg$correlation_filter_threshold %||% 0.95
        # remove columns with too many NAs or zero var
        na_frac_allowed <- cfg$preproc$max_na_fraction %||% 0.2
        valid_cols <- sapply(as.data.frame(Xtr_scaled), function(col) mean(is.na(col)) <= na_frac_allowed & sd(col, na.rm=TRUE) > 0)
        if (sum(valid_cols) == 0) next
        Xtr_for_corr <- as.data.frame(Xtr_scaled[, valid_cols, drop = FALSE])
        sel_after_corr <- correlation_feature_filter(Xtr_for_corr, threshold = corr_threshold)
        # subset both Xtr_scaled and corresponding Xte_scaled to selected columns
        # if correlation filter returns empty, fall back to all valid cols
        if (length(sel_after_corr) == 0) sel_after_corr <- colnames(Xtr_for_corr)
        Xtr_scaled_filt <- Xtr_scaled[, sel_after_corr, drop = FALSE]
        Xte_scaled_filt <- Xte_scaled[, sel_after_corr, drop = FALSE]
        if (ncol(Xtr_scaled_filt) == 0) next

        # Ranking: primary is LASSO; optional Ranger-based ranking can be triggered by config
        fea_list <- cfg$fea_list %||% c(1,3,5,7,9,10,12,15,20,25,30)
        fea_list <- sort(unique(as.integer(fea_list[fea_list > 0])))

        # LASSO ranking on inner-train filtered features
        cvfit <- cv.glmnet(x = Xtr_scaled_filt, y = y_train, family = "binomial", alpha = cfg$lasso$alpha %||% 1, nfolds = cfg$lasso$nfolds %||% 5, type.measure = "auc", parallel = FALSE)
        lasso_coefs <- as.numeric(coef(cvfit, s = "lambda.min")[-1]); names(lasso_coefs) <- colnames(Xtr_scaled_filt)
        lasso_rank <- names(sort(abs(lasso_coefs), decreasing = TRUE))

        # Ranger ranking: only if configured for ranking (importance_for_ranking), otherwise skip
        importance_for_ranking <- cfg$ranger$importance_for_ranking %||% cfg$ranger$importance %||% "none"
        if (importance_for_ranking %in% c("impurity","permutation")) {
          y_train_factor <- factor(y_train, levels = c(0,1))
          rfit_rank <- ranger::ranger(x = as.data.frame(Xtr_scaled_filt), y = y_train_factor, probability = TRUE, importance = importance_for_ranking, num.trees = cfg$ranger$num.trees %||% 500)
          ranger_imp <- rfit_rank$variable.importance
          ranger_rank <- names(sort(ranger_imp, decreasing = TRUE))
        } else {
          ranger_rank <- character(0)
        }

        # For each k: select features and compute AUCs for both methods
        for (k in fea_list) {
          k_use <- min(k, ncol(Xtr_scaled_filt))
          if (k_use <= 0) next

          # LASSO top-k selection
          auc_lasso_k <- NA_real_
          sel_lasso <- intersect(lasso_rank, colnames(Xtr_scaled_filt))[seq_len(k_use)]
          if (length(sel_lasso) > 0) {
            lres <- fit_predict_lasso_safe(Xtr_scaled_filt[, sel_lasso, drop = FALSE], y_train, Xte_scaled_filt[, sel_lasso, drop = FALSE], y_test, nfolds = cfg$lasso$nfolds %||% 5, alpha = cfg$lasso$alpha %||% 1)
            auc_lasso_k <- lres$auc
          }

          # Ranger top-k selection: use ranger_rank if available; final ranger model trained with importance = "none"
          auc_ranger_k <- NA_real_
          sel_ranger <- if (length(ranger_rank) > 0) intersect(ranger_rank, colnames(Xtr_scaled_filt))[seq_len(k_use)] else character(0)
          # if user explicitly requested Ranger ranking = none, but still wants Ranger prediction for top-k,
          # we fallback to use LASSO ranking for selecting features for Ranger prediction (keeps consistent selection)
          if (length(sel_ranger) == 0 && cfg$ranger$use_lasso_for_ranger_selection %||% TRUE) {
            sel_ranger <- sel_lasso
          }
          if (length(sel_ranger) > 0) {
            Xtr_r_k <- as.data.frame(Xtr_scaled_filt[, sel_ranger, drop = FALSE])
            Xte_r_k <- as.data.frame(Xte_scaled_filt[, sel_ranger, drop = FALSE])
            # final ranger trained with importance = "none" to avoid using importance for prediction
            rfit_k <- ranger::ranger(x = Xtr_r_k, y = factor(y_train, levels = c(0,1)), probability = TRUE, importance = "none", num.trees = cfg$ranger$num.trees %||% 500)
            pred_obj_k <- predict(rfit_k, data = Xte_r_k)
            preds_ranger_k <- get_ranger_prob1(pred_obj_k, rfit_k)
            auc_ranger_k <- tryCatch(as.numeric(pROC::auc(y_test, preds_ranger_k)), error = function(e) NA_real_)
          }

          # append both methods explicitly
          new_rows <- data.table(cohort = cohort_name, mode = mode, outer = o, fold = fold_id, method = "LASSO", n_features = k_use, auc = auc_lasso_k)
          new_rows <- rbind(new_rows, data.table(cohort = cohort_name, mode = mode, outer = o, fold = fold_id, method = "Ranger", n_features = k_use, auc = auc_ranger_k))
          all_progression_dt <- rbind(all_progression_dt, new_rows)

          # periodic flush to disk to bound memory if many folds
          if (nrow(all_progression_dt) >= flush_threshold) {
            flush_part_idx <- flush_part_idx + 1
            flush_progression_to_disk(all_progression_dt, outdir, flush_part_idx)
            all_progression_dt <- all_progression_dt[0]
          }
        } # end k loop

      } # end inner folds
    } # end outer fold loop
  } # end mode loop
} # end cohort loop

# flush remaining progression if any and collect into single table
if (nrow(all_progression_dt) > 0) {
  flush_part_idx <- flush_part_idx + 1
  flush_progression_to_disk(all_progression_dt, outdir, flush_part_idx)
  # read back all partials to build unified table
  partials <- list.files(file.path(outdir, "metrics"), pattern = "^progression_partial_.*\\.csv$", full.names = TRUE)
  all_parts <- lapply(partials, fread)
  all_progression_dt_full <- rbindlist(all_parts, use.names = TRUE, fill = TRUE)
  # save unified
  fwrite(all_progression_dt_full, file.path(outdir, "metrics", "all_progression_full.csv"))
} else {
  # if nothing flushed but maybe earlier flushes exist, read them
  partials <- list.files(file.path(outdir, "metrics"), pattern = "^progression_partial_.*\\.csv$", full.names = TRUE)
  if (length(partials) == 0) stop("No progression results collected.")
  all_parts <- lapply(partials, fread)
  all_progression_dt_full <- rbindlist(all_parts, use.names = TRUE, fill = TRUE)
}

# ---------------- Cross-cohort evaluation (train on full cohort -> predict on other cohort(s)) ----------------
if (isTRUE(cfg$cross_eval %||% FALSE)) {
  cohort_names <- names(cohort_data_list)
  for (i in seq_along(cohort_names)) {
    for (j in seq_along(cohort_names)) {
      if (i == j) next
      train_name <- cohort_names[i]; test_name <- cohort_names[j]
      cat(sprintf("\nCross-eval: train %s -> test %s\n", train_name, test_name))
      train_dat <- cohort_data_list[[train_name]]$tsv; test_dat <- cohort_data_list[[test_name]]$tsv
      meta_keep <- c("ID", cfg$case_column %||% "CASE_CONTROL", "SEX"); meta_keep <- meta_keep[meta_keep %in% colnames(train_dat)]
      feats_train <- setdiff(colnames(train_dat), meta_keep); feats_test <- setdiff(colnames(test_dat), meta_keep)
      common_feats <- intersect(feats_train, feats_test)
      if (length(common_feats) == 0) { cat("No common features for cross-eval\n"); next }
      X_train_full <- as.matrix(train_dat[, common_feats, drop = FALSE]); y_train_full <- as.numeric(train_dat[[cfg$case_column %||% "CASE_CONTROL"]])
      X_test_full <- as.matrix(test_dat[, common_feats, drop = FALSE]); y_test_full <- as.numeric(test_dat[[cfg$case_column %||% "CASE_CONTROL"]])
      # scale using train
      train_center <- colMeans(X_train_full, na.rm = TRUE); train_scale <- apply(X_train_full, 2, sd, na.rm = TRUE); train_scale[train_scale == 0] <- 1
      Xtr_scaled <- scale(X_train_full, center = train_center, scale = train_scale); Xte_scaled <- scale(X_test_full, center = train_center, scale = train_scale)
      # correlation filter on train_full
      valid_cols <- sapply(as.data.frame(Xtr_scaled), function(col) sd(col, na.rm = TRUE) > 0)
      Xtr_corr <- as.data.frame(Xtr_scaled[, valid_cols, drop = FALSE])
      sel_after_corr <- correlation_feature_filter(Xtr_corr, threshold = cfg$preproc$correlation_filter_threshold %||% 0.95)
      if (length(sel_after_corr) == 0) sel_after_corr <- colnames(Xtr_corr)
      Xtr_scaled_filt <- Xtr_scaled[, sel_after_corr, drop = FALSE]; Xte_scaled_filt <- Xte_scaled[, sel_after_corr, drop = FALSE]
      # ranking on train_full
      cvf <- cv.glmnet(x = Xtr_scaled_filt, y = y_train_full, family = "binomial", alpha = cfg$lasso$alpha %||% 1, nfolds = cfg$lasso$nfolds %||% 5, type.measure = "auc")
      lasso_coefs_full <- as.numeric(coef(cvf, s = "lambda.min")[-1]); names(lasso_coefs_full) <- colnames(Xtr_scaled_filt)
      lasso_rank_full <- names(sort(abs(lasso_coefs_full), decreasing = TRUE))
      importance_for_ranking <- cfg$ranger$importance_for_ranking %||% cfg$ranger$importance %||% "none"
      if (importance_for_ranking %in% c("impurity","permutation")) {
        rfit_rank_full <- ranger::ranger(x = as.data.frame(Xtr_scaled_filt), y = factor(y_train_full, levels = c(0,1)), probability = TRUE, importance = importance_for_ranking, num.trees = cfg$ranger$num.trees %||% 500)
        ranger_rank_full <- names(sort(rfit_rank_full$variable.importance, decreasing = TRUE))
      } else ranger_rank_full <- character(0)
      fea_list <- cfg$fea_list %||% c(1,3,5,7,9,10,12,15,20,25,30)
      for (k in fea_list) {
        k_use <- min(k, ncol(Xtr_scaled_filt)); if (k_use <= 0) next
        sel_l <- intersect(lasso_rank_full, colnames(Xtr_scaled_filt))[seq_len(k_use)]
        sel_r <- if (length(ranger_rank_full) > 0) intersect(ranger_rank_full, colnames(Xtr_scaled_filt))[seq_len(k_use)] else character(0)
        if (length(sel_r) == 0 && cfg$ranger$use_lasso_for_ranger_selection %||% TRUE) sel_r <- sel_l
        # LASSO predict on test
        auc_l <- NA_real_
        if (length(sel_l) > 0) {
          lres <- fit_predict_lasso_safe(Xtr_scaled_filt[, sel_l, drop = FALSE], y_train_full, Xte_scaled_filt[, sel_l, drop = FALSE], y_test_full, nfolds = cfg$lasso$nfolds %||% 5, alpha = cfg$lasso$alpha %||% 1)
          auc_l <- lres$auc
        }
        # Ranger predict on test
        auc_r <- NA_real_
        if (length(sel_r) > 0) {
          rfit_k <- ranger::ranger(x = as.data.frame(Xtr_scaled_filt[, sel_r, drop = FALSE]), y = factor(y_train_full, levels = c(0,1)), probability = TRUE, importance = "none", num.trees = cfg$ranger$num.trees %||% 500)
          pred_obj <- predict(rfit_k, data = as.data.frame(Xte_scaled_filt[, sel_r, drop = FALSE]))
          preds_r <- get_ranger_prob1(pred_obj, rfit_k)
          auc_r <- tryCatch(as.numeric(pROC::auc(y_test_full, preds_r)), error = function(e) NA_real_)
        }
        all_progression_dt_full <- data.table(cohort = paste0(train_name,"->",test_name), mode = "cross", outer = NA_integer_, fold = NA_integer_, method = "LASSO", n_features = k_use, auc = auc_l)
        all_progression_dt_full <- rbind(all_progression_dt_full, data.table(cohort = paste0(train_name,"->",test_name), mode = "cross", outer = NA_integer_, fold = NA_integer_, method = "Ranger", n_features = k_use, auc = auc_r))
        # write these rows to disk and later combine with per-fold ones
        f_cross <- file.path(outdir, "metrics", paste0("cross_progression_", train_name, "_to_", test_name, ".csv"))
        if (!file.exists(f_cross)) fwrite(all_progression_dt_full, f_cross) else fwrite(all_progression_dt_full, f_cross, append = TRUE)
      }
    }
  }
}

# ----------------------- Aggregation and final summary plots -----------------------

# read back all progressions: per-fold partials + cross files
partials <- list.files(file.path(outdir, "metrics"), pattern = "^progression_partial_.*\\.csv$", full.names = TRUE)
perfolds <- if (length(partials) > 0) rbindlist(lapply(partials, fread), use.names = TRUE, fill = TRUE) else data.table()
cross_files <- list.files(file.path(outdir, "metrics"), pattern = "^cross_progression_.*\\.csv$", full.names = TRUE)
cross_dt <- if (length(cross_files) > 0) rbindlist(lapply(cross_files, fread), use.names = TRUE, fill = TRUE) else data.table()
# also include all_progression_dt_full if exists
if (exists("all_progression_dt_full")) {
  all_dt <- rbindlist(list(perfolds, all_progression_dt_full, cross_dt), use.names = TRUE, fill = TRUE)
} else {
  all_dt <- rbindlist(list(perfolds, cross_dt), use.names = TRUE, fill = TRUE)
}

if (nrow(all_dt) == 0) stop("No progression results available for aggregation.")

# aggregate across outer folds for each cohort/mode/method/n_features
agg_dt <- all_dt[, .(
  mean_auc = mean(auc, na.rm = TRUE),
  sd_auc = sd(auc, na.rm = TRUE),
  n = sum(!is.na(auc))
), by = .(cohort, mode, method, n_features)]
agg_dt[, ci_lower := mean_auc - 1.96 * (sd_auc / sqrt(pmax(n,1)))]
agg_dt[, ci_upper := mean_auc + 1.96 * (sd_auc / sqrt(pmax(n,1)))]

# save aggregated table
fwrite(agg_dt, file.path(outdir, "metrics", "agg_auc_vs_features.csv"))

# produce final plots per cohort/mode
settings <- unique(agg_dt[, .(cohort, mode)])
final_plots <- list()
for (i in seq_len(nrow(settings))) {
  c <- settings$cohort[i]; m <- settings$mode[i]
  sub <- agg_dt[cohort == c & mode == m]
  if (nrow(sub) == 0) next
  p <- ggplot(sub, aes(x = n_features, y = mean_auc, color = method, group = method)) +
    geom_line(size = 1) + geom_point(size = 2) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = method), alpha = 0.15, color = NA) +
    scale_x_continuous(breaks = sort(unique(sub$n_features))) +
    labs(title = paste0("AUC vs #features  ", c, " (", m, ")"), x = "Number of top features (k)", y = "Mean AUC �95% CI") +
    theme_minimal() + theme(legend.position = "bottom")
  ggsave(file.path(outdir, "figures", paste0("auc_vs_features_summary_", gsub("[/:>]","_",c), "_", m, ".pdf")), p, width = 8, height = 6)
  final_plots[[length(final_plots) + 1]] <- p
}

# multi-page PDF if requested
if (isTRUE(cfg$io$combine_plots_pdf %||% FALSE) && length(final_plots) > 0) {
  pdf(file.path(outdir, "figures", "all_auc_vs_features_summary.pdf"), width = 8, height = 6)
  for (g in final_plots) print(g)
  dev.off()
}

stopCluster(cl)
cat("Pipeline finished. Results in:", outdir, "\n")
