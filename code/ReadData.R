 #!/usr/bin/env Rscript
# ReadData.R - robustes Einlesen + Vorbereitung für Cohorts (bo, bg, ...)
# - liest config.yaml
# - findet und liest .fam, .eigenvec, .tsv
# - matched IDs, erstellt CASE_CONTROL falls nötig (aus .fam)
# - hängt SEX + PCs an
# - optionale Residualisierung der PRS auf PCs+SEX
# - einfache PopStrat (PC outlier detection)
#
# Benötigt Pakete: data.table, yaml, stats, broom (optional)
# Install: install.packages(c("data.table","yaml","broom"))

library(data.table)
library(yaml)
library(stats)   # lm
library(broom)   # optional für tidy() (nice diagnostics)

# -----------------------------
# Hilfsfunktionen
# -----------------------------

# sichere Datei-Suche: erstes passendes file oder stop mit Fehlermeldung
find_first_file <- function(dir, pattern) {
  files <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) stop(sprintf("Keine Dateien gefunden in %s mit Muster %s", dir, pattern))
  return(files[1])
}

# erkennt ID-Spalte im TSV anhand Kandidatenliste, sonst index/erste Spalte
detect_tsv_id_col <- function(dt, id_candidates) {
  for (c in id_candidates) {
    if (c %in% colnames(dt)) return(c)
  }
  # fallback: falls index text ist -> reset index, name it ID
  if (!is.null(rownames(dt)) && any(rownames(dt) != 1:nrow(dt))) {
    dt2 <- copy(dt)
    dt2[, ID__index := rownames(dt)]
    setcolorder(dt2, c("ID__index", setdiff(names(dt2), "ID__index")))
    setnames(dt2, "ID__index", "ID")
    return("ID")
  }
  # sonst erste Spalte
  return(colnames(dt)[1])
}

# map PLINK .fam PHENO zu 0/1 (control=0, case=1), accepts plink codes 1/2 or 0/1
map_fam_pheno_to_binary <- function(pheno_vec) {
  # pheno_vec numeric
  pheno_vec[pheno_vec == -9] <- NA   # missing
  # if values are 1/2 -> map 1->0 control, 2->1 case
  if (all(pheno_vec %in% c(1,2,NA))) {
    mapped <- ifelse(pheno_vec == 2, 1, ifelse(pheno_vec == 1, 0, NA))
    return(mapped)
  }
  # if 0/1 already
  if (all(pheno_vec %in% c(0,1,NA))) return(pheno_vec)
  # else try to coerce: treat >1 as case
  return(ifelse(is.na(pheno_vec), NA, ifelse(pheno_vec > 1, 1, 0)))
}

# residualize matrix df_prs (data.table/data.frame) on covariates (data.frame)
residualize_prs <- function(df_prs, covariates) {
  # df_prs: data.frame/data.table with PRS columns
  # covariates: data.frame with same number of rows, covariate columns (including intercept if desired)
  # returns data.frame of residuals, same colnames
  if (nrow(df_prs) != nrow(covariates)) stop("residualize_prs: len mismatch")
  resid_mat <- as.data.frame(matrix(NA, nrow = nrow(df_prs), ncol = ncol(df_prs)))
  colnames(resid_mat) <- colnames(df_prs)
  # build model matrix once
  X <- model.matrix(~ ., data = covariates)
  XtX_inv <- tryCatch(solve(crossprod(X)), error = function(e) NULL)
  # faster OLS using lm for each column (robust enough)
  for (j in seq_len(ncol(df_prs))) {
    y <- df_prs[[j]]
    if (all(is.na(y))) {
      resid_mat[[j]] <- NA
      next
    }
    # check constant column
    if (var(y, na.rm = TRUE) == 0) {
      resid_mat[[j]] <- 0  # constant -> zero residual
      next
    }
    fit <- lm(y ~ ., data = covariates)
    resid_mat[[j]] <- resid(fit)
  }
  return(resid_mat)
}

# detect PC outliers: returns logical mask of non-outlier rows
pc_outlier_mask <- function(pcs_df, sd_threshold = 6, pcs_to_use = c("PC1","PC2")) {
  # pcs_df: data.frame with PC1..PCn
  # outlier if abs(z) > sd_threshold on any of the pcs_to_use
  stopifnot(all(pcs_to_use %in% colnames(pcs_df)))
  mask <- rep(TRUE, nrow(pcs_df))
  for (pc in pcs_to_use) {
    z <- scale(pcs_df[[pc]])
    mask <- mask & (abs(z) <= sd_threshold)
  }
  return(mask)
}

# -----------------------------
# Hauptteil: lese Config und verarbeite Cohorts
# -----------------------------
cfg_file <- "config.yaml"
if (!file.exists(cfg_file)) stop("config.yaml nicht gefunden im Arbeitsverzeichnis")

cfg <- yaml.load_file(cfg_file)

out_dir <- cfg$output_dir %||% "outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

id_candidates <- cfg$id_candidates %||% c("ID","IID","Sample_ID")
case_col <- cfg$case_column %||% "CASE_CONTROL"
pcs_prefix <- cfg$pcs_prefix %||% "PC"
use_resid <- isTRUE(cfg$use_residualize)
resid_covs <- cfg$resid_covariates %||% c("SEX", paste0("PC",1:10))
outlier_sd <- cfg$outlier_sd %||% 6
n_pcs_use <- cfg$n_pcs_use %||% 10

for (cohort in cfg$cohorts) {
  name <- cohort$name
  cat(sprintf("\n--- Processing cohort: %s ---\n", name))

  fam_dir <- cohort$fam_eig_dir
  tsv_dir <- cohort$tsv_dir
  tsv_pattern <- cohort$tsv_pattern %||% "\\.tsv$"

  # find files
  fam_file <- find_first_file(fam_dir, "\\.fam$")
  eig_file <- find_first_file(fam_dir, "\\.eigenvec$")
  tsv_file <- find_first_file(tsv_dir, tsv_pattern)

  cat("FAM:", fam_file, "\n")
  cat("EIGENVEC:", eig_file, "\n")
  cat("TSV:", tsv_file, "\n")

  # read files
  fam <- fread(fam_file, header = FALSE, data.table = FALSE)   # PLINK .fam (6 cols)
  eig <- fread(eig_file, header = FALSE, data.table = FALSE)   # often first two cols are FID/IID then PCs
  tsv <- fread(tsv_file, header = TRUE, data.table = FALSE)    # multiPRS table

  # detect id col in tsv
  tsv_id_col <- detect_tsv_id_col(tsv, id_candidates)
  if (tsv_id_col != "ID") {
    if ("ID" %in% colnames(tsv)) {
      # drop existing ID to avoid duplicates, then rename chosen
      tsv$ID <- NULL
    }
    setnames(tsv, old = tsv_id_col, new = "ID")
  }
  tsv$ID <- as.character(tsv$ID)

  # ensure fam/eig ids are character
  fam[[2]] <- as.character(fam[[2]])
  eig[[2]] <- as.character(eig[[2]])

  # CASE_CONTROL detection: if not in tsv, take from fam (col 6)
  if (!(case_col %in% colnames(tsv))) {
    if (ncol(fam) >= 6) {
      fam_pheno <- data.frame(ID = fam[[2]], PHENO = as.numeric(fam[[6]]), stringsAsFactors = FALSE)
      fam_pheno$CASE_FROM_FAM <- map_fam_pheno_to_binary(fam_pheno$PHENO)
      # merge into tsv on ID (safe merge)
      if ("CASE_FROM_FAM" %in% colnames(tsv)) tsv$CASE_FROM_FAM <- NULL
      tsv <- merge(tsv, fam_pheno[, c("ID","CASE_FROM_FAM")], by = "ID", all.x = TRUE)
      tsv[[case_col]] <- tsv$CASE_FROM_FAM
      tsv$CASE_FROM_FAM <- NULL
      cat(sprintf("Info: '%s' nicht in TSV; Phenotype aus .fam gemerged in '%s'.\n", case_col, case_col))
    } else {
      stop(sprintf("Kein Fall/Control in TSV und .fam hat keine phenotype (col6). Für %s Abbruch.", name))
    }
  }

  # compute common IDs
  common_ids <- Reduce(intersect, list(fam[[2]], eig[[2]], tsv$ID))
  cat(sprintf("Common IDs count for %s: %d\n", name, length(common_ids)))
  if (length(common_ids) == 0) stop("Keine gemeinsamen IDs!")

  # filter data to common ids and keep order stable
  fam_common <- fam[fam[[2]] %in% common_ids, , drop = FALSE]
  eig_common <- eig[eig[[2]] %in% common_ids, , drop = FALSE]
  tsv_common <- tsv[tsv$ID %in% common_ids, , drop = FALSE]

  # reorder by ID to be consistent (optional)
  # we'll create an index order based on tsv_common$ID
  order_ids <- tsv_common$ID
  fam_common <- fam_common[match(order_ids, fam_common[[2]]), , drop = FALSE]
  eig_common <- eig_common[match(order_ids, eig_common[[2]]), , drop = FALSE]
  tsv_common <- tsv_common[match(order_ids, tsv_common$ID), , drop = FALSE]

  # rename/eig PCP columns: assume eig columns are (FID, IID, PC1, PC2, ...), sometimes first col is dropped -> handle both
  # If eig has >2 cols, build PC names
  n_eig_cols <- ncol(eig_common)
  if (n_eig_cols >= 3) {
    pc_count <- min(n_eig_cols - 2, n_pcs_use)
    pc_names <- paste0(pcs_prefix, 1:pc_count)
    colnames(eig_common)[1:(2 + pc_count)] <- c("FID","ID", pc_names)
    pcs_df <- as.data.frame(eig_common[, c("ID", pc_names)])
  } else {
    # not expected but fallback: no PCs
    pcs_df <- data.frame(ID = eig_common[[2]])
    pc_names <- character(0)
    cat("Warnung: keine PCs in eigenvec gefunden.\n")
  }

  # SEX: from fam col 5
  if (ncol(fam_common) >= 5) {
    sex_vec <- fam_common[[5]]
    # PLINK: sex 1=male,2=female; keep as-is or map as needed
    tsv_common$SEX <- sex_vec
  } else {
    cat("Warnung: keine SEX-Spalte in .fam (V5) gefunden.\n")
    tsv_common$SEX <- NA
  }

  # add PCs to tsv_common by ID merge (IDs aligned but merge to be safe)
  if (nrow(pcs_df) > 0) {
    tsv_common <- merge(tsv_common, pcs_df, by = "ID", all.x = TRUE, sort = FALSE)
  }

  # optional: detect PRS columns automatically: numeric columns excluding ID,CASE_CONTROL,SEX,PCs
  exclude_cols <- c("ID", case_col, "SEX", pc_names)
  prs_cols <- setdiff(colnames(tsv_common)[sapply(tsv_common, is.numeric)], exclude_cols)
  # remove columns with single unique value
  prs_cols <- prs_cols[sapply(prs_cols, function(cn) length(unique(tsv_common[[cn]][!is.na(tsv_common[[cn]])])) > 1)]
  cat(sprintf("Detected %d numeric PRS-like columns for %s (sample): %s\n", length(prs_cols), name, paste(head(prs_cols,10), collapse=", ")))

  # Optional residualization: regress each PRS on covariates (e.g., PCs + SEX) and keep residuals
  if (use_resid && length(prs_cols) > 0) {
    covs_here <- resid_covs[resid_covs %in% colnames(tsv_common)]
    if (length(covs_here) == 0) {
      cat("Warnung: keine Residual-Kovariaten gefunden. Residualisierung übersprungen.\n")
    } else {
      cat(sprintf("Residualisiere %d PRS-Spalten auf Kovariaten: %s\n", length(prs_cols), paste(covs_here, collapse=", ")))
      cov_df <- tsv_common[, covs_here, drop = FALSE]
      # convert SEX to numeric if factor-like
      if ("SEX" %in% colnames(cov_df)) cov_df$SEX <- as.numeric(as.character(cov_df$SEX))
      # call residualize_prs (returns a data.frame)
      resid_mat <- residualize_prs(tsv_common[, prs_cols, drop = FALSE], cov_df)
      # replace in tsv_common with residuals, suffix to indicate residualized
      for (j in seq_along(prs_cols)) {
        tsv_common[[prs_cols[j]]] <- resid_mat[[j]]
      }
      cat("Residualisierung abgeschlossen.\n")
    }
  } else {
    cat("Residuum-Schritt deaktiviert oder keine PRS-Spalten vorhanden.\n")
  }

  # population stratification: detect outliers on PC1/PC2 (z-score threshold)
  if (all(c("PC1","PC2") %in% colnames(tsv_common))) {
    mask_keep <- pc_outlier_mask(tsv_common[, c("PC1","PC2"), drop = FALSE], sd_threshold = outlier_sd, pcs_to_use = c("PC1","PC2"))
    removed <- sum(!mask_keep)
    cat(sprintf("PC outlier detection: removed %d samples (threshold %g SD)\n", removed, outlier_sd))
    if (removed > 0) {
      tsv_common <- tsv_common[mask_keep, , drop = FALSE]
      fam_common <- fam_common[mask_keep, , drop = FALSE]
      eig_common <- eig_common[mask_keep, , drop = FALSE]
    }
  } else {
    cat("Warnung: PC1/PC2 nicht vorhanden – Outlier-Removal übersprungen.\n")
  }

  # final: sanity checks
  tsv_common <- tsv_common[order(tsv_common$ID), ]
  stopifnot(!any(duplicated(tsv_common$ID)))
  cat(sprintf("Final dataset for %s: %d samples x %d cols\n", name, nrow(tsv_common), ncol(tsv_common)))
# nach dem Erzeugen von tsv_common, aber vor write.table()
# ermittle zuerst deine PRS-Spalten (wie eh schon gemacht)
  exclude_cols <- c("ID", case_col, "SEX", pc_names)
  prs_cols     <- setdiff(colnames(tsv_common), exclude_cols)

# gewünschte neue Reihenfolge
  new_order <- c("ID", case_col, "SEX", pc_names, prs_cols)

# setze die Spalten in dieser Reihenfolge
  tsv_common <- tsv_common[, new_order]

  # write prepared file
  out_file <- file.path(out_dir, sprintf("%s_prepared.tsv", name))
  write.table(tsv_common, out_file, sep = "\t", row.names = FALSE, quote = FALSE)
  cat(sprintf("Saved prepared file: %s\n", out_file))
}

