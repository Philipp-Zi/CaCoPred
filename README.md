# CaCoPred
Its a bit wild west in the beginning, but should evolve into a well structured git for CC-Pred in MultiPRS files

## Input Data

The project requires three types of input files:

- **.fam**
  Contains sample IDs and phenotypes (case/control status).

- **.eigen**
  Contains eigenvectors from principal component analysis (PCA), used to correct for population stratification.

- **.tsv**
  Contains polygenic risk scores (PRS) for each individual.

### File Organization
All input files should be placed in the `data/` directory.
Make sure filenames are consistent with the configuration specified in `config/default.yaml` or `config/local.yaml`.

# cclogreg – Case-Control Prediction with one PRS and one Cohord

## Step 1: Single-PRS Case-Control Logistic Regression

This script (`cclogreg.R`) performs a logistic regression using a single polygenic risk score (PRS) to predict case-control status.

### Input
- **PRS file** (`*.best` or `*.score` from PRSice)
  Must contain columns: `FID`, `IID`, `PRS`
- **FAM file**
  Must contain columns: `FID`, `IID`, `PHENO`
  - `PHENO = 1` → control
  - `PHENO = 2` → case

### Output
- `output/cc_logreg_results.tsv` → summary table with sample size, odds ratio, confidence interval, and AUC
- `output/cc_logreg_ROC.png` → ROC curve plot

### Example Run
```bash
Rscript cclogreg.R
