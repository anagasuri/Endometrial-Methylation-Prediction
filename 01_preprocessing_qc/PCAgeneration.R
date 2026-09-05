############################################################
# PRE vs POST SmartSVA PCA plots
# Color by binary cycle phase:
#   1) Proliferative
#   2) Secretory
#
# This script does NOT run SmartSVA again.
# It:
#   - rebuilds PRE matrix from raw data
#   - loads POST corrected matrix from RDS
#   - runs PCA on each
#   - collapses cycle phase to Proliferative vs Secretory
#   - saves 2 new PCA plots
#   - prints PE/SE sample counts
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
M_RDS         <- "M_batch1_batch2.rds"
KEEP_IDS_TXT <- "R01_Study_IDs_to_include_220622.txt"
ANNOT_CSV    <- "SH-Data Annotation_07.07.20_1.csv"
POST_RDS     <- "M_mat_SmartSVA_corrected.rds"

# -----------------------------
# Output directory
# -----------------------------
OUT_DIR <- "/mnt/efs/home/ubuntu/endo-me_data/PCA_case_cycle_pre_post/PE_SE_only_plots"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 1) Load raw methylation data
# -----------------------------
M <- readRDS(M_RDS)

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

# -----------------------------
# 5) Extract barcodes
# -----------------------------
barcodes <- trimws(annot_sub$`Epic_Complete Bar code`)

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

# -----------------------------
# 7) Convert to numeric PRE matrix
# -----------------------------
cpg_ids <- M_sub[[1]]
M_mat <- as.matrix(M_sub[, -1, drop = FALSE])
rownames(M_mat) <- cpg_ids

rm(M, M_sub)
gc()

# -----------------------------
# 8) Align covariates to PRE matrix columns
# -----------------------------
covariates <- as.data.frame(annot_sub)
rownames(covariates) <- barcodes
covariates <- covariates[colnames(M_mat), , drop = FALSE]

stopifnot(identical(rownames(covariates), colnames(M_mat)))

# -----------------------------
# 9) Collapse cycle phase to Proliferative vs Secretory
# -----------------------------
covariates$`Cycle phase for Analysis` <- trimws(
  as.character(covariates$`Cycle phase for Analysis`)
)

covariates$Cycle_binary <- ifelse(
  covariates$`Cycle phase for Analysis` == "PE",
  "Proliferative",
  "Secretory"
)

covariates$Cycle_binary <- factor(
  covariates$Cycle_binary,
  levels = c("Proliferative", "Secretory")
)

cycle_binary_colors <- c(
  "Proliferative" = "blue",
  "Secretory" = "orange"
)

# -----------------------------
# 10) Run PRE PCA
# -----------------------------
pca_pre <- prcomp(t(M_mat), center = FALSE, scale. = FALSE)
write.csv(pca_pre$x, file.path(OUT_DIR, "pca_scores_pre_PE_SE_only.csv"))

# -----------------------------
# 11) Load POST corrected matrix
# -----------------------------
M_corrected <- readRDS(POST_RDS)

M_corrected <- M_corrected[, colnames(M_mat), drop = FALSE]

stopifnot(identical(colnames(M_corrected), colnames(M_mat)))
stopifnot(identical(colnames(M_corrected), rownames(covariates)))

# -----------------------------
# 12) Run POST PCA
# -----------------------------
pca_post <- prcomp(t(M_corrected), center = FALSE, scale. = FALSE)
write.csv(pca_post$x, file.path(OUT_DIR, "pca_scores_post_PE_SE_only.csv"))

# -----------------------------
# 13) Plot helper
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
  plot(
    pca_obj$x[, 1],
    pca_obj$x[, 2],
    pch = 16,
    col = cols,
    xlab = "PC1",
    ylab = "PC2",
    main = main_txt
  )
  legend(
    "topright",
    legend = levels(labs),
    col = legend_cols,
    pch = 16,
    cex = 0.7
  )
  dev.off()
}

# -----------------------------
# 14) Save PRE/POST binary cycle phase PCA plots
# -----------------------------
plot_pc(
  pca_pre,
  covariates$Cycle_binary,
  "Pre-SmartSVA PCA: Proliferative vs Secretory",
  file.path(OUT_DIR, "PCA_pre_proliferative_vs_secretory.png"),
  colors = cycle_binary_colors
)

plot_pc(
  pca_post,
  covariates$Cycle_binary,
  "Post-SmartSVA PCA: Proliferative vs Secretory",
  file.path(OUT_DIR, "PCA_post_proliferative_vs_secretory.png"),
  colors = cycle_binary_colors
)

# -----------------------------
# 15) Quick checks and sample counts
# -----------------------------
cat("\nWorking directory:\n")
print(getwd())

cat("\nPRE matrix dimensions:\n")
print(dim(M_mat))

cat("\nPOST matrix dimensions:\n")
print(dim(M_corrected))

cat("\nOriginal cycle phase counts:\n")
print(table(covariates$`Cycle phase for Analysis`, useNA = "ifany"))

cat("\nCollapsed Proliferative vs Secretory counts:\n")
print(table(covariates$Cycle_binary, useNA = "ifany"))

cat("\nOutput files written to:\n")
print(OUT_DIR)

cat("\nFiles in output directory:\n")
print(list.files(OUT_DIR))

# -----------------------------
# -----------------------------
# -----------------------------

# -----------------------------
# Count PE vs Secretory samples only
# Exclude Menstrual, NA, unknown, or unclear labels
# -----------------------------

phase <- trimws(as.character(covariates$`Cycle phase for Analysis`))

cycle_group <- rep(NA_character_, length(phase))

# Proliferative
cycle_group[phase == "PE"] <- "Proliferative"

# Secretory
cycle_group[phase %in% c("ESE", "MSE", "LSE", "SE")] <- "Secretory"

# Keep only clearly labeled Proliferative or Secretory samples
keep_cycle <- !is.na(cycle_group)

cycle_counts <- table(
  factor(cycle_group[keep_cycle], levels = c("Proliferative", "Secretory"))
)

excluded_counts <- table(phase[!keep_cycle], useNA = "ifany")

cat("\nOriginal cycle phase counts:\n")
print(table(phase, useNA = "ifany"))

cat("\nIncluded PE vs Secretory counts:\n")
print(cycle_counts)

cat("\nNumber of included samples:\n")
print(sum(cycle_counts))

cat("\nExcluded samples/counts:\n")
print(excluded_counts)


