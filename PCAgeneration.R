############################################################
# PRE vs POST SmartSVA PCA plots
# Color by:
#   1) Endometriosis (Yes/No)
#   2) Cycle phase for Analysis
#
# This script does NOT run SmartSVA again.
# It:
#   - rebuilds PRE matrix from raw data
#   - loads POST corrected matrix from RDS
#   - runs PCA on each
#   - saves 4 plots
############################################################

# -----------------------------
# Libraries
# -----------------------------
suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------------
# Input files
# -----------------------------
M_RDS        <- "M_batch1_batch2.rds"
KEEP_IDS_TXT <- "R01_Study_IDs_to_include_220622.txt"
ANNOT_CSV    <- "SH-Data Annotation_07.07.20_1.csv"
POST_RDS     <- "M_mat_SmartSVA_corrected.rds"

# -----------------------------
# Output directory
# -----------------------------
OUT_DIR <- "/home/ubuntu/endo-me_data/PCA_case_cycle_pre_post"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 1) Load raw methylation data
# -----------------------------
M <- readRDS(M_RDS)   # CpGs in col1, sample barcodes in other cols

# -----------------------------
# 2) Read keep IDs
# -----------------------------
keep_ids <- readLines(KEEP_IDS_TXT)
keep_ids <- trimws(keep_ids)
keep_ids <- keep_ids[keep_ids != "Study.ID"]

# -----------------------------
# 3) Read annotation CSV
# -----------------------------
annot <- fread(ANNOT_CSV, data.table = FALSE)
colnames(annot) <- trimws(colnames(annot))

# -----------------------------
# 4) Subset annotation to keep_ids
# -----------------------------
annot$`Study ID` <- trimws(annot$`Study ID`)
annot_sub <- annot[annot$`Study ID` %in% keep_ids, ]
nrow(annot_sub)

# -----------------------------
# 5) Extract barcodes
# -----------------------------
barcodes <- trimws(annot_sub$`Epic_Complete Bar code`)
length(barcodes)
head(barcodes)

# -----------------------------
# 6) Subset methylation matrix to those barcodes
# -----------------------------
idx <- match(barcodes, colnames(M))

if (any(is.na(idx))) {
  missing_barcodes <- barcodes[is.na(idx)]
  stop("These barcodes are in annotation but missing from M colnames:\n",
       paste(missing_barcodes, collapse = "\n"))
}

M_sub <- M[, c(1, idx), drop = FALSE]
dim(M_sub)

# -----------------------------
# 7) Convert to numeric PRE matrix
# -----------------------------
cpg_ids <- M_sub[[1]]
M_mat <- as.matrix(M_sub[, -1, drop = FALSE])
rownames(M_mat) <- cpg_ids

dim(M_mat)

# free memory
rm(M, M_sub)
gc()

# -----------------------------
# 8) Align covariates to PRE matrix columns
# -----------------------------
covariates <- as.data.frame(annot_sub)
rownames(covariates) <- barcodes
covariates <- covariates[colnames(M_mat), , drop = FALSE]

stopifnot(identical(rownames(covariates), colnames(M_mat)))

# make factors
covariates$`Endometriosis (Yes/No)` <- factor(covariates$`Endometriosis (Yes/No)`)
covariates$`Cycle phase for Analysis` <- factor(covariates$`Cycle phase for Analysis`)

# -----------------------------
# 8.5) Define colors
# -----------------------------
case_colors <- c(
  "No"  = "black",
  "Yes" = "red"
)

# -----------------------------
# 9) Run PRE PCA
# -----------------------------
pca_pre <- prcomp(t(M_mat), center = FALSE, scale. = FALSE)
write.csv(pca_pre$x, file.path(OUT_DIR, "pca_scores_pre.csv"))

# -----------------------------
# 10) Load POST corrected matrix
# -----------------------------
M_corrected <- readRDS(POST_RDS)

# force same sample order
M_corrected <- M_corrected[, colnames(M_mat), drop = FALSE]

stopifnot(identical(colnames(M_corrected), colnames(M_mat)))
stopifnot(identical(colnames(M_corrected), rownames(covariates)))

# -----------------------------
# 11) Run POST PCA
# -----------------------------
pca_post <- prcomp(t(M_corrected), center = FALSE, scale. = FALSE)
write.csv(pca_post$x, file.path(OUT_DIR, "pca_scores_post.csv"))

# -----------------------------
# 12) Plot helper
# -----------------------------
plot_pc <- function(pca_obj, labels, main_txt, out_file, colors = NULL) {
  labs <- as.factor(labels)
  
  if (is.null(colors)) {
    cols <- rainbow(length(levels(labs)))[labs]
    legend_cols <- rainbow(length(levels(labs)))
  } else {
    cols <- colors[as.character(labs)]
    legend_cols <- colors[levels(labs)]
  }
  
  png(out_file, width = 2000, height = 1600, res = 300)
  plot(pca_obj$x[,1], pca_obj$x[,2],
       pch = 16,
       col = cols,
       xlab = "PC1",
       ylab = "PC2",
       main = main_txt)
  legend("topright",
         legend = levels(labs),
         col = legend_cols,
         pch = 16,
         cex = 0.7)
  dev.off()
}

# -----------------------------
# 13) Save PRE plots
# -----------------------------
plot_pc(
  pca_pre,
  covariates$`Endometriosis (Yes/No)`,
  "Pre-SmartSVA PCA: Endometriosis (Yes/No)",
  file.path(OUT_DIR, "PCA_pre_case_control.png"),
  colors = case_colors
)

# plot_pc(
#   pca_pre,
#   covariates$`Cycle phase for Analysis`,
#   "Pre-SmartSVA PCA: Cycle phase for Analysis",
#   file.path(OUT_DIR, "PCA_pre_cycle_phase.png")
# )

# -----------------------------
# 14) Save POST plots
# -----------------------------
plot_pc(
  pca_post,
  covariates$`Endometriosis (Yes/No)`,
  "Post-SmartSVA PCA: Endometriosis (Yes/No)",
  file.path(OUT_DIR, "PCA_post_case_control.png"),
  colors = case_colors
)

# plot_pc(
#   pca_post,
#   covariates$`Cycle phase for Analysis`,
#   "Post-SmartSVA PCA: Cycle phase for Analysis",
#   file.path(OUT_DIR, "PCA_post_cycle_phase.png")
# )

# -----------------------------
# 15) Quick checks
# -----------------------------
print(getwd())
print(dim(M_mat))
print(dim(M_corrected))
print(list.files(OUT_DIR))