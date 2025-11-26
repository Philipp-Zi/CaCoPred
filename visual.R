#!/usr/bin/env Rscript

# 1. Libraries
library(jsonlite)
library(data.table)
library(ggplot2)
library(reshape2)
library(stringr)

# 2. Pfade zu den Resultaten
res_file   <- "outputs/all_cohorts_cv_summary.json"
topk_file  <- "outputs/top_features_log.json"

# 3. Daten laden
summary_list <- fromJSON(res_file, simplifyDataFrame = FALSE)
topk_log     <- fromJSON(topk_file)

# 4. Metriken in Data.Frame umformen
# Für jedes cohort und jedes k extrahieren wir Medianwerte
plot_df <- rbindlist(lapply(names(summary_list), function(cohort) {
  df <- summary_list[[cohort]]
  dt <- rbindlist(lapply(names(df), function(k) {
    d0 <- df[[k]]
    data.table(cohort = cohort,
               k       = as.integer(k),
               AUC     = d0$AUC_median,
               OR      = d0$OR_median,
               COR     = d0$COR_median,
               R2      = d0$R2_median)
  }))
  dt
}))

# 5. Lineplot: Median AUC, OR, R2 vs k
for(metric in c("AUC","OR","COR","R2")) {
  ggplot(plot_df, aes(x = k, y = get(metric), color = cohort)) +
    geom_line() + geom_point() +
    scale_x_continuous(breaks = unique(plot_df$k)) +
    ggtitle(sprintf("Median %s vs Top-k-Features", metric)) +
    xlab("k (Anzahl Features)") + ylab(metric) +
    theme_minimal() -> p
  print(p)
}

# 6. Top-k-Feature-Stabilität
# Tabelle mit Feature-Auswahl-Frequenzen pro k
# key in topk_log: "<cohort>_k_<k>_fold_<f>_outer_<o>"
entries <- lapply(names(topk_log), function(key) {
  parts <- str_split(key, "_")[[1]]
  cohort <- parts[1]
  k      <- as.integer(parts[3])
  fv     <- topk_log[[key]]
  data.table(cohort = cohort, k = k, feature = fv)
})
feat_dt <- rbindlist(entries)

# Zähle wie oft jedes Feature in topk (über Folds × Outers) auftaucht
freq_dt <- feat_dt[, .(count = .N),
                   by = .(cohort, k, feature)]

# 7. Top-10 Features nach Auswahlhäufigkeit pro k
top_feats <- freq_dt[order(-count)][, head(.SD, 10), by = .(cohort, k)]

# 8. Barplot: Auswahl-häufigkeit der 10 häufigsten Features (für ein cohort & k)
for(c in unique(top_feats$cohort)) {
  for(K in unique(top_feats$k)) {
    sub <- top_feats[cohort == c & k == K]
    ggplot(sub, aes(x = reorder(feature, count), y = count)) +
      geom_col(fill = "steelblue") +
      coord_flip() +
      ggtitle(sprintf("Top 10 Features nach Häufigkeit\ncohort=%s, k=%d", c, K)) +
      ylab("Auswahl-Frequenz") + xlab("Feature") +
      theme_minimal() -> b
    print(b)
  }
}

# 9. Heatmap: Feature-Auswahl über k-Werte
heat_dt <- freq_dt[, .(freq = count), by = .(feature, k)]
heat_wide <- dcast(heat_dt, feature ~ k, value.var = "freq", fill = 0)
mat <- as.matrix(heat_wide[,-1])
rownames(mat) <- heat_wide$feature
library(pheatmap)
pheatmap(mat,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         main = "Feature-Stabilität über k")

# 10. Speichern (optional)
# ggsave("auc_vs_k.png", p)
# saveRDS(plot_df, "plot_df.RDS")
