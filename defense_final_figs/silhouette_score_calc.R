############################################################
# Silhouette scores from existing PRE/POST PCA scores
# Full dataset, no PCA rerun, no SmartSVA rerun
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(cluster)
  library(ggplot2)
})

############################################################
# FILE PATHS
############################################################

PCA_PRE  <- "/home/ubuntu/endo-me_data/PCA_case_cycle_pre_post/pca_scores_pre.csv"
PCA_POST <- "/home/ubuntu/endo-me_data/PCA_case_cycle_pre_post/pca_scores_post.csv"

ANNOT_CSV    <- "/mnt/efs/home/ubuntu/SH-Data Annotation_07.07.20_1.csv"
KEEP_IDS_TXT <- "/mnt/efs/home/ubuntu/R01_Study_IDs_to_include_220622.txt"

OUT_DIR <- "/mnt/efs/home/ubuntu/defense_final_figs/"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

############################################################
# LOAD METADATA
############################################################

keep_ids <- readLines(KEEP_IDS_TXT)
keep_ids <- trimws(keep_ids)
keep_ids <- keep_ids[keep_ids != "Study.ID"]

annot <- fread(ANNOT_CSV, data.table = FALSE)
colnames(annot) <- trimws(colnames(annot))

annot$`Study ID` <- trimws(annot$`Study ID`)
annot_sub <- annot[annot$`Study ID` %in% keep_ids, ]

annot_sub$sample_id <- trimws(annot_sub$`Epic_Complete Bar code`)

############################################################
# LOAD PCA SCORES
############################################################

pca_pre <- read.csv(PCA_PRE, check.names = FALSE)
pca_post <- read.csv(PCA_POST, check.names = FALSE)

# First column is sample ID
colnames(pca_pre)[1] <- "sample_id"
colnames(pca_post)[1] <- "sample_id"

############################################################
# FUNCTION TO CALCULATE SILHOUETTE
############################################################

calc_silhouette <- function(pca_df, meta_df, label_col, n_pcs = 10) {
  
  df <- pca_df %>%
    inner_join(meta_df, by = "sample_id")
  
  pc_cols <- grep("^PC", colnames(df), value = TRUE)
  pc_cols <- pc_cols[1:min(n_pcs, length(pc_cols))]
  
  df <- df %>%
    filter(!is.na(.data[[label_col]]))
  
  labels <- as.factor(df[[label_col]])
  
  # Need at least 2 groups
  if (length(levels(labels)) < 2) {
    return(NA)
  }
  
  X <- df[, pc_cols, drop = FALSE]
  
  sil <- silhouette(as.numeric(labels), dist(X))
  mean(sil[, 3])
}

############################################################
# CALCULATE SCORES
############################################################

results <- data.frame(
  label = c(
    "Institute for Analysis",
    "Batch",
    "Endometriosis (Yes/No)",
    "Cycle phase for Analysis"
  ),
  pre = c(
    calc_silhouette(pca_pre, annot_sub, "Institute for Analysis"),
    calc_silhouette(pca_pre, annot_sub, "Batch"),
    calc_silhouette(pca_pre, annot_sub, "Endometriosis (Yes/No)"),
    calc_silhouette(pca_pre, annot_sub, "Cycle phase for Analysis")
  ),
  post = c(
    calc_silhouette(pca_post, annot_sub, "Institute for Analysis"),
    calc_silhouette(pca_post, annot_sub, "Batch"),
    calc_silhouette(pca_post, annot_sub, "Endometriosis (Yes/No)"),
    calc_silhouette(pca_post, annot_sub, "Cycle phase for Analysis")
  )
)

results$change <- results$post - results$pre

print(results)

write.csv(
  results,
  file.path(OUT_DIR, "silhouette_scores_pre_post.csv"),
  row.names = FALSE
)

############################################################
# PLOT RESULTS
############################################################

plot_df <- results %>%
  tidyr::pivot_longer(
    cols = c(pre, post),
    names_to = "SmartSVA_status",
    values_to = "silhouette_score"
  )

plot_df$SmartSVA_status <- factor(
  plot_df$SmartSVA_status,
  levels = c("pre", "post"),
  labels = c("Pre-SmartSVA", "Post-SmartSVA")
)

p <- ggplot(plot_df, aes(x = label, y = silhouette_score, fill = SmartSVA_status)) +
  geom_col(position = position_dodge(width = 0.8)) +
  theme_bw(base_size = 12) +
  labs(
    title = "Silhouette Scores Before and After SmartSVA",
    x = "",
    y = "Mean silhouette score",
    fill = ""
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

ggsave(
  file.path(OUT_DIR, "silhouette_scores_pre_post.png"),
  p,
  width = 9,
  height = 5,
  dpi = 300
)

cat("Saved results to:\n", OUT_DIR, "\n")