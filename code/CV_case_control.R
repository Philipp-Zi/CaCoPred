# code/CV_case_control.R

#–– Libraries ––#
library(yaml)
library(data.table)
library(caret)
library(glmnet)
library(ranger)
library(pROC)
library(jsonlite)

# null-coalescing operator
`%||%` <- function(a, b) if (!is.null(a)) a else b

#–– Helper‐Funktionen ––#

sex_stratified_mask <- function(group_vec, sex_vec) {
  stopifnot(length(group_vec) == length(sex_vec))
  df       <- data.frame(g = factor(group_vec),
                         s = factor(sex_vec),
                         idx = seq_along(group_vec))
  keep_idx <- integer(0)
  for (gval in levels(df$g)) {
    sub    <- df[df$g == gval, ]
    counts <- table(sub$s)
    if (length(counts) == 1) {
      keep_idx <- c(keep_idx, sub$idx)
    } else {
      minc <- min(counts)
      for (sexv in names(counts)) {
        chosen <- sample(sub$idx[sub$s == sexv], size = minc, replace = FALSE)
        keep_idx <- c(keep_idx, chosen)
      }
    }
  }
  mask <- rep(FALSE, length(group_vec))
  mask[keep_idx] <- TRUE
  mask
}

select_best_threshold_per_gwas <- function(train_df, prs_cols, label_vec,
  gwas_extract = function(col) tail(strsplit(col, "_")[[1]], 1)
) {
  gwas_ids   <- sapply(prs_cols, gwas_extract)
  unique_gwas <- unique(gwas_ids)
  chosen     <- character(0)
  for (g in unique_gwas) {
    cols <- prs_cols[gwas_ids == g]
    cors <- vapply(cols, function(cn) {
      x <- train_df[[cn]]
      if (length(unique(na.omit(x))) <= 1) return(0)
      r <- suppressWarnings(cor(x, label_vec, use = "complete.obs"))
      if (is.na(r)) 0 else abs(r)
    }, FUN.VALUE = numeric(1))
    chosen <- c(chosen, cols[which.max(cors)])
  }
  chosen
}

nagelkerke_r2 <- function(y, p) {
  if (is.matrix(p) && ncol(p) > 1) p <- p[, 2]
  y <- as.numeric(y)
  if (length(unique(y)) == 1) return(NA)
  m0 <- glm(y ~ 1, family = binomial)
  m1 <- glm(y ~ p, family = binomial)
  ll0 <- as.numeric(logLik(m0))
  ll1 <- as.numeric(logLik(m1))
  n   <- length(y)
  r2_cs  <- 1 - exp((2 * (ll0 - ll1)) / n)
  r2_max <- 1 - exp((2 * ll0) / n)
  if (r2_max == 0) NA else r2_cs / r2_max
}

majority_vote_labels <- function(pred_list) {
  mat <- do.call(rbind, pred_list)
  apply(mat, 2, function(col) {
    ux <- sort(table(col), decreasing = TRUE)
    as.numeric(names(ux)[1])
  })
}

median_proba <- function(proba_list) {
  mat <- do.call(rbind, proba_list)
  apply(mat, 2, median, na.rm = TRUE)
}

#–– Config & Setup ––#
cfg       <- yaml::read_yaml("config.yaml")
outdir    <- cfg$output_dir      %||% "outputs"
case_col  <- cfg$case_column     %||% "CASE_CONTROL"
sex_col   <- cfg$sex_column      %||% "SEX"
pcs_pref  <- cfg$pc_prefix       %||% "PC"
fea_list  <- cfg$fea_list        %||% c(1,3,5,10,20)
innerloop <- as.integer(cfg$innerloop   %||% 5)
nun_iter  <- as.integer(cfg$nun_iter    %||% 3)
n_splits  <- as.integer(cfg$n_splits    %||% 5)
seed      <- as.integer(cfg$random_seed %||% 42)
n_estim   <- as.integer(cfg$n_estimators %||% 500)
mode      <- cfg$prs_mode        %||% "best_threshold_per_gwas"

set.seed(seed)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

results_all <- list()

#–– Main Loop über Kohorten ––#
for (cohort in cfg$cohorts) {
  name      <- cohort$name
  prep_file <- file.path(outdir, paste0(name, "_prepared.tsv"))
  if (!file.exists(prep_file)) stop("Prepared file not found: ", prep_file)
  df <- fread(prep_file, data.table = FALSE)

  # PRS-Spalten ermitteln
  exclude_cols <- c("ID", case_col, sex_col,
                    grep(sprintf("^%s", pcs_pref), colnames(df), value = TRUE))
  prs_cols <- setdiff(colnames(df)[sapply(df, is.numeric)], exclude_cols)
  if (length(prs_cols) == 0) stop("No PRS columns in ", name)

  cat(sprintf("\n=== Cohort %s ===\nDetected %d PRS columns\n",
              name, length(prs_cols)))

  # Speicher für Ergebnisse pro k
  cohort_res <- lapply(fea_list, function(k) {
    list(OR = numeric(), COR = numeric(), R2 = numeric(), AUC = numeric())
  })
  names(cohort_res) <- as.character(fea_list)

  # Outer Repeats
  for (outer in seq_len(nun_iter)) {
    cat(sprintf(" Outer iter %d/%d\n", outer, nun_iter))
    folds <- createFolds(df[[case_col]], k = n_splits)

    for (fold_idx in seq_along(folds)) {
      cat(sprintf("  Fold %d/%d\n", fold_idx, length(folds)))
      test_idx      <- folds[[fold_idx]]
      train_df_full <- df[-test_idx, , drop = FALSE]
      test_df_full  <- df[ test_idx, , drop = FALSE]

      # Per Feature-Size k
      for (k in fea_list) {
        inner_labels <- list()
        inner_probs  <- list()

        for (inner in seq_len(innerloop)) {
          mask     <- sex_stratified_mask(train_df_full[[case_col]],
                                          train_df_full[[sex_col]])
          train_df <- train_df_full[mask, , drop = FALSE]

          # Y definieren
          y_raw   <- train_df[[case_col]]
          stopifnot(all(na.omit(y_raw) %in% c(0,1)))
          y_train <- droplevels(factor(y_raw, levels = c(0,1)))
          y_test  <- as.numeric(test_df_full[[case_col]])

          # Default-Features: alle PRS + PCs
          feat_cols <- c(prs_cols,
                         grep(sprintf("^%s", pcs_pref),
                              colnames(df),
                              value = TRUE))

          if (length(levels(y_train)) < 2) {
            # Single-Class → konstante Vorhersage
            only_cls   <- as.numeric(levels(y_train))[1]
            pred_proba <- rep(only_cls, length(y_test))
            pred_label <- rep(only_cls, length(y_test))

          } else if (mode == "best_threshold_per_gwas") {
            # 1) Pro-GWAS Top-Threshold wählen
            chosen <- select_best_threshold_per_gwas(
              train_df, prs_cols, as.numeric(y_train)
            )

            # 2) Punkt-biseriale Korrelationen absolut berechnen
            cor_vals <- vapply(
              chosen,
              FUN = function(cn) {
                x <- train_df[[cn]]
                if (length(unique(na.omit(x))) <= 1) return(0)
                r <- suppressWarnings(cor(x, as.numeric(y_train), use = "complete.obs"))
                if (is.na(r)) 0 else abs(r)
              },
              FUN.VALUE = numeric(1)
            )

            # 3) Top-k bestimmen
            ord       <- order(cor_vals, decreasing = TRUE)
            topk      <- chosen[ord][seq_len(min(k, length(ord)))]
            feat_cols <- c(
              topk,
              grep(sprintf("^%s", pcs_pref), colnames(df), value = TRUE)
            )

            # 4) Debug: wie viele Features, was fehlt
            missing_feats <- setdiff(feat_cols, colnames(test_df_full))
            cat(sprintf("   >> Fold %d, k=%d → %d Features; missing in test: %d\n",
                        fold_idx, k, length(feat_cols), length(missing_feats)))
            if (length(missing_feats) > 0) {
              cat("      fehlende Features: ", paste(missing_feats, collapse = ", "), "\n")
            }

                        # 5) TrainX bauen (robust)
            trainX <- as.data.frame(train_df[, feat_cols, drop = FALSE])
            # ensure column names are syntactically valid and unique (consistent with test)
            colnames(trainX) <- make.names(colnames(trainX), unique = TRUE)

            # build testX defensively: initialize zeros for all feat_cols and then fill with real values
            testX <- as.data.frame(matrix(0,
                                          nrow = nrow(test_df_full),
                                          ncol = length(feat_cols),
                                          dimnames = list(NULL, make.names(feat_cols, unique = TRUE))))
            # fill where available
            available <- intersect(colnames(test_df_full), feat_cols)
            if (length(available) > 0) {
              # use make.names mapping to match column names used in trainX/testX
              mapped_names <- setNames(make.names(available, unique = TRUE), available)
              for (orig in available) {
                safe_name <- make.names(orig, unique = TRUE)
                # assign, but ensure length matches
                testX[[safe_name]] <- as.numeric(test_df_full[[orig]])
              }
            }

            # Reorder testX to match exact trainX column order
            # If some train columns are missing in testX (shouldn't after initialization), add them as 0
            missing_cols <- setdiff(colnames(trainX), colnames(testX))
            if (length(missing_cols) > 0) {
              # add missing cols filled with 0
              for (mc in missing_cols) testX[[mc]] <- 0
            }
            # ensure same column order
            testX <- testX[, colnames(trainX), drop = FALSE]

            # Coerce all to numeric (train may contain integer/factor)
            trainX[] <- lapply(trainX, function(x) as.numeric(as.character(x)))
            testX[]  <- lapply(testX,  function(x) as.numeric(as.character(x)))

            # final sanity check
            missing_in_test <- setdiff(colnames(trainX), colnames(testX))
            if (length(missing_in_test) > 0) {
              cat("CRITICAL: Still missing in test (skipping this inner run):\n")
              print(missing_in_test)
              # push NA predictions (so the fold contributes NA metrics) and continue
              pred_proba <- rep(NA_real_, nrow(test_df_full))
              pred_label <- rep(NA_real_, nrow(test_df_full))
            } else {
              # Fit model (Y as factor with levels 0/1)
              # ensure response vector matches rows of trainX
              Y <- train_df[[case_col]]
              # If Y is numeric 0/1, convert to factor with labels "0","1" (ranger makes columns "0"/"1")
              if (is.numeric(Y)) Y <- factor(Y, levels = c(0,1))
              if (is.character(Y)) Y <- factor(Y) # hope levels are 0/1

              rf <- tryCatch(
                ranger(dependent.variable.name = "Y",
                       data = data.frame(Y = Y, trainX),
                       probability = TRUE,
                       num.trees = n_estim,
                       num.threads = 1),
                error = function(e) e
              )
              if (inherits(rf, "error")) {
                cat("ranger training error, skipping this inner run:\n"); print(rf)
                pred_proba <- rep(NA_real_, nrow(test_df_full))
                pred_label <- rep(NA_real_, nrow(test_df_full))
              } else {
                pr <- tryCatch(predict(rf, data = testX)$predictions, error = function(e) e)
                if (inherits(pr, "error")) {
                  cat("ranger predict error, skipping this inner run:\n"); print(pr)
                  pred_proba <- rep(NA_real_, nrow(test_df_full))
                  pred_label <- rep(NA_real_, nrow(test_df_full))
                } else {
                  # pr might be matrix with colnames "0","1" or "predictions" column; pick col "1"
                  if ("1" %in% colnames(pr)) {
                    pred_proba <- pr[, "1"]
                  } else if (ncol(pr) >= 2) {
                    pred_proba <- pr[, ncol(pr)]
                  } else {
                    # fallback: single column probabilities or votes
                    pred_proba <- as.numeric(pr[,1])
                  }
                  pred_label <- as.numeric(pred_proba >= 0.5)
                }
              }
            }


          } else if (mode == "glmnet_all_thresholds") {
            # glmnet branch unverändert
            feat_cols <- feat_cols  # all thresholds = keine Selektion
            x_train   <- as.matrix(train_df[, feat_cols, drop = FALSE])
            x_test    <- as.matrix(test_df_full[, feat_cols, drop = FALSE])
            cvfit     <- cv.glmnet(x_train,
                                   as.numeric(y_train),
                                   family       = "binomial",
                                   alpha        = 0.5,
                                   nfolds       = 5,
                                   type.measure = "auc")
            pred_proba <- predict(cvfit,
                                  newx = x_test,
                                  type = "response",
                                  s    = "lambda.min")[, 1]
            pred_label <- as.numeric(pred_proba >= 0.5)

          } else {
            stop("Unknown mode: ", mode)
          }

          inner_labels[[length(inner_labels) + 1]] <- pred_label
          inner_probs [[length(inner_probs)  + 1]] <- pred_proba
        }  # inner

        # –– Aggregation & Metriken ––
        agg_l   <- majority_vote_labels(inner_labels)
        agg_p   <- median_proba(inner_probs)
        ct      <- table(y_test, agg_l)
        or_val  <- if (all(dim(ct) == c(2,2))) tryCatch(fisher.test(ct)$estimate,
                                                        error = function(e) NA) else NA
        cor_val <- if (length(unique(agg_l))>1 && length(unique(y_test))>1)
                     tryCatch(cor(as.numeric(agg_l), y_test),
                              error = function(e) NA) else NA
        r2_val  <- tryCatch(nagelkerke_r2(y_test, agg_p), error = function(e) NA)
        auc_val <- tryCatch(auc(roc(y_test, agg_p)), error = function(e) NA)

        cohort_res[[as.character(k)]]$OR  <- c(cohort_res[[as.character(k)]]$OR,  or_val)
        cohort_res[[as.character(k)]]$COR <- c(cohort_res[[as.character(k)]]$COR, cor_val)
        cohort_res[[as.character(k)]]$R2  <- c(cohort_res[[as.character(k)]]$R2,  r2_val)
        cohort_res[[as.character(k)]]$AUC <- c(cohort_res[[as.character(k)]]$AUC, auc_val)
      }  # k
    }    # fold
  }      # outer

  #–– Cohort-Summary ––#
  summary_list <- lapply(names(cohort_res), function(k) {
    arr_or  <- unlist(cohort_res[[k]]$OR)
    arr_cor <- unlist(cohort_res[[k]]$COR)
    arr_r2  <- unlist(cohort_res[[k]]$R2)
    arr_auc <- unlist(cohort_res[[k]]$AUC)
    list(
      OR_median  = median(arr_or,  na.rm = TRUE),
      COR_median = median(arr_cor, na.rm = TRUE),
      R2_median  = median(arr_r2,  na.rm = TRUE),
      AUC_median = median(arr_auc, na.rm = TRUE),
      OR_all     = arr_or,
      COR_all    = arr_cor,
      R2_all     = arr_r2,
      AUC_all    = arr_auc
    )
  })
  names(summary_list) <- names(cohort_res)

  saveRDS(summary_list,
          file = file.path(outdir, paste0(name, "_cv_summary.RDS")))
  write_json(summary_list,
             path      = file.path(outdir, paste0(name, "_cv_summary.json")),
             pretty    = TRUE, auto_unbox = TRUE)
  results_all[[name]] <- summary_list
}

#–– Final Save ––#
saveRDS(results_all,
        file = file.path(outdir, "all_cohorts_cv_summary.RDS"))
write_json(results_all,
           path      = file.path(outdir, "all_cohorts_cv_summary.json"),
           pretty    = TRUE, auto_unbox = TRUE)

cat("\nAll done.\n")
