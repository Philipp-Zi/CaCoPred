#!/usr/bin/env Rscript

library(data.table)
library(pROC)

# ---- INPUT ----
# Replace these with the paths to your PRS and FAM files
prs_file <- "input/cohort_prs.all_score"   # PRSice output file
fam_file <- "input/cohort_data.fam"        # FAM file with phenotype

# ---- Load Data ----
prs <- fread(prs_file, header = TRUE)
fam <- fread(fam_file, header = FALSE)

# Standard FAM columns: FID IID PID MID SEX PHENO
colnames(fam) <- c("FID","IID","PID","MID","SEX","PHENO")

# Merge PRS + phenotype
dat <- merge(prs, fam[, .(FID, IID, PHENO)], by = c("FID","IID"))

# Recode phenotype: 1 = control, 2 = case
dat$Status <- ifelse(dat$PHENO == 2, 1,
                     ifelse(dat$PHENO == 1, 0, NA))

# ---- Logistic Regression ----
# Select the desired PRS column, e.g. Pt_0.05
model <- glm(Status ~ Pt_0.05, data = dat, family = binomial())

# Odds Ratio and Confidence Interval
or <- exp(coef(model)["Pt_0.05"])
ci <- exp(confint(model, "Pt_0.05"))

# Predictions and ROC
preds <- predict(model, type = "response")
roc_obj <- roc(dat$Status, preds)

# ---- Console Output ----
cat("Sample size:", nrow(dat), "\n")
cat("Odds Ratio:", round(or, 3),
    "95% CI:", round(ci[1], 3), "-", round(ci[2], 3), "\n")
cat("AUC:", round(auc(roc_obj), 3), "\n")

# ---- Save Results ----
out_dir <- "output"
if (!dir.exists(out_dir)) dir.create(out_dir)

results <- data.frame(
  Cohort = "example_cohort",
  N = nrow(dat),
  OR = round(or, 3),
  CI_low = round(ci[1], 3),
  CI_high = round(ci[2], 3),
  AUC = round(auc(roc_obj), 3)
)

fwrite(results, file.path(out_dir, "cc_logreg_results.tsv"), sep = "\t")

png(file.path(out_dir, "cc_logreg_ROC.png"))
plot(roc_obj, main = "ROC Curve for PRS @ 0.05")
dev.off()
