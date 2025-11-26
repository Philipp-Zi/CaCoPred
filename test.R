#!/usr/bin/env Rscript
library(jsonlite)
library(stringr)

outdir <- "outputs"
cv_file <- file.path(outdir, "all_cohorts_cv_summary.json")
feat_file <- file.path(outdir, "top_features_log.json")

cat("WD:", getwd(), "\n")
cat("Sizes: CV=", file.size(cv_file), " FEAT=", file.size(feat_file), "\n\n")

cv_raw <- fromJSON(cv_file, simplifyVector = FALSE)
feat_raw <- fromJSON(feat_file, simplifyVector = FALSE)

cat("cv_raw: class:", class(cv_raw), " length:", length(cv_raw), "\n")
cat("First 20 keys cv_raw:\n"); print(head(names(cv_raw),20))
cat("\nfeat_raw: class:", class(feat_raw), " length:", length(feat_raw), "\n")
cat("First 20 keys feat_raw:\n"); print(head(names(feat_raw),20))

cat("\nStructure cv_raw[[1]]:\n"); str(cv_raw[[1]], max.level = 2)
cat("\nStructure feat_raw[[1]]:\n"); str(feat_raw[[1]], max.level = 2)

# Quick scan for AUC-containing keys
extract_auc_try <- function(x) {
  if (is.list(x)) {
    poss <- c("auc","AUC","mean_auc","meanAUC","cv_auc","roc_auc","mean")
    for (p in poss) if (!is.null(x[[p]])) return(TRUE)
    nums <- suppressWarnings(as.numeric(unlist(x)))
    return(length(nums[!is.na(nums)])>0)
  }
  if (is.numeric(x)) return(TRUE)
  return(FALSE)
}
has_auc <- sapply(cv_raw, extract_auc_try)
cat("\nKeys with AUC-like content:", sum(has_auc), "of", length(has_auc), "\n")
print(head(names(cv_raw)[which(has_auc)], 20))

# Quick scan for features being character-like
get_feat_try <- function(x) {
  if (is.character(x)) return(TRUE)
  if (is.list(x)) {
    poss <- c("top","features","selected","top_features")
    for (p in poss) if (!is.null(x[[p]])) return(TRUE)
    # any nested character values?
    allch <- unlist(x)
    return(any(sapply(allch, is.character)))
  }
  return(FALSE)
}
has_feat <- sapply(feat_raw, get_feat_try)
cat("\nKeys with feature-like content:", sum(has_feat), "of", length(has_feat), "\n")
print(head(names(feat_raw)[which(has_feat)], 20))

