#–– Libraries ––#
library(yaml)
library(data.table)

`%||%` <- function(a, b) if (!is.null(a)) a else b

#–– Helpers ––#

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

#–– Main ––#

cfg_file <- "config.yaml"
if (!file.exists(cfg_file)) stop("config.yaml not found in working directory")
cfg <- yaml::yaml.load_file(cfg_file)

out_dir       <- cfg$output_dir      %||% "outputs"
id_candidates <- cfg$id_candidates   %||% c("ID","IID","Sample_ID")
case_col      <- cfg$case_column     %||% "CASE_CONTROL"
pcs_prefix    <- cfg$pc_prefix       %||% "PC"
n_pcs_use     <- cfg$n_pcs_use       %||% 10
outlier_sd    <- cfg$outlier_sd      %||% 6

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

for (cohort in cfg$cohorts) {
  print(cfg$cohorts)
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

  tsv_id_col <- detect_tsv_id_col(tsv, id_candidates)
  if (tsv_id_col != "ID") {
    if ("ID" %in% colnames(tsv)) tsv$ID <- NULL
    setnames(tsv, tsv_id_col, "ID")
  }
  tsv$ID <- as.character(tsv$ID)

  # Add CASE_CONTROL if missing (from fam col 6)
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

  # SEX from fam col 5 (PLINK: 1=male, 2=female), then filter invalid
  if (!("SEX" %in% colnames(tsv))) {
    if (ncol(fam) >= 5) {
      sex_vec <- as.numeric(fam[[5]])
      tsv <- merge(tsv, data.frame(ID = fam[[2]], SEX = sex_vec), by = "ID", all.x = TRUE)
    } else {
      tsv$SEX <- NA_real_
      warning("SEX not found in .fam; set to NA.")
    }
  }
  # Ensure numeric and filter rows with SEX not in {1,2}
  tsv$SEX <- as.numeric(tsv$SEX)
  valid_sex_mask <- !is.na(tsv$SEX) & tsv$SEX %in% c(1, 2)
  removed_sex <- sum(!valid_sex_mask)
  if (removed_sex > 0) {
    cat(sprintf("SEX filtering: removed %d rows with SEX not in {1,2}.\n", removed_sex))
    tsv <- tsv[valid_sex_mask, , drop = FALSE]
  }

  # PCs
  pcs_df <- extract_pcs(eig, pcs_prefix = pcs_prefix, n_pcs_use = n_pcs_use)
  tsv <- merge(tsv, pcs_df, by = "ID", all.x = TRUE)

  # Common IDs and ordering
  common_ids <- Reduce(intersect, list(fam[[2]], eig[[2]], tsv$ID))
  cat(sprintf("Common IDs: %d\n", length(common_ids)))
  if (length(common_ids) == 0) stop("No common IDs!")
  tsv <- tsv[tsv$ID %in% common_ids, , drop = FALSE]
  tsv <- tsv[order(tsv$ID), , drop = FALSE]

  # Optional PC outlier removal (PC1/PC2)
  pc1 <- paste0(pcs_prefix, 1); pc2 <- paste0(pcs_prefix, 2)
  if (all(c(pc1, pc2) %in% colnames(tsv))) {
    mask <- pc_outlier_mask(tsv[, c(pc1, pc2)], sd_threshold = outlier_sd)
    removed <- sum(!mask)
    if (removed > 0) {
      cat(sprintf("PC outlier removal: %d samples removed (>%g SD).\n", removed, outlier_sd))
      tsv <- tsv[mask, , drop = FALSE]
    }
  } else {
    cat("PC1/PC2 not present; skipping outlier removal.\n")
  }

  stopifnot(!any(duplicated(tsv$ID)))
  cat(sprintf("Final dataset for %s: %d samples x %d cols\n", name, nrow(tsv), ncol(tsv)))
#–– Feature-Filter und Spaltenreihenfolge ––#

# Hole Präfix aus Config oder setze Default
feature_prefix <- cfg$feature_prefix %||% "Pt_1"

# Alle gewünschten Features mit dem Präfix
feature_cols <- grep(paste0("^", feature_prefix), colnames(tsv), value = TRUE)

# Alle anderen Spalten, die NICHT Feature-Spalten sind
meta_cols <- colnames(tsv)[!grepl("^Pt_", colnames(tsv))]

# Neue Reihenfolge: erst Meta, dann gefilterte Features
tsv <- tsv[, c(meta_cols, feature_cols), drop = FALSE]

cat("Kept", length(feature_cols), "features with prefix", feature_prefix, "and moved them to the end\n")

  # No residualization here (to avoid leakage). Provide raw PRS + SEX + PCs + CASE_CONTROL.
  out_file <- file.path(out_dir, paste0(name, "_prepared.tsv"))
  fwrite(tsv, out_file, sep = "\t", quote = FALSE, na = "NA")
  cat("Written: ", out_file, "\n")
}

cat("\nAll cohorts prepared with strict SEX filtering and without leakage-prone residualization.\n")
