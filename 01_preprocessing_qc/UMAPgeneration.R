############################################################
# UMAP from existing PCA scores
#
# Creates UMAP plots from PRE and POST PCA scores.
#
# Outputs:
#   - Case/control UMAP plots
#   - PE vs Secretory UMAP plots
#
# Cycle phase grouping:
#   PE  -> Proliferative
#   ESE, MSE, LSE, SE -> Secretory
#   Menstrual/NA/unclear labels are excluded from PE vs Secretory plots
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(uwot)
})

# -----------------------------
# Input files
# -----------------------------
PCA_PRE  <- "/mnt/efs/home/ubuntu/endo-me_data/PCA_case_cycle_pre_post/pca_scores_pre.csv"
PCA_POST <- "/mnt/efs/home/ubuntu/endo-me_data/PCA_case_cycle_pre_post/pca_scores_post.csv"

ANNOT_CSV    <- "SH-Data Annotation_07.07.20_1.csv"
KEEP_IDS_TXT <- "R01_Study_IDs_to_include_220622.txt"

# -----------------------------
# Output directories
# -----------------------------
OUT_DIR <- "/mnt/efs/home/ubuntu/endo-me_data/UMAP_from_PCA/PE_SE_only_plots_umaps"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

OUT_DIR_BINARY <- "/mnt/efs/home/ubuntu/endo-me_data/UMAP_from_PCA/PE_SE_only_plots_umaps"
dir.create(OUT_DIR_BINARY, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Parameters
# -----------------------------
N_PCS <- 30

# -----------------------------
# Load PCA scores
# -----------------------------
pca_pre  <- fread(PCA_PRE, data.table = FALSE)
pca_post <- fread(PCA_POST, data.table = FALSE)

rownames(pca_pre)  <- trimws(pca_pre[, 1])
rownames(pca_post) <- trimws(pca_post[, 1])

pca_pre  <- pca_pre[, -1, drop = FALSE]
pca_post <- pca_post[, -1, drop = FALSE]

rownames(pca_pre)  <- trimws(rownames(pca_pre))
rownames(pca_post) <- trimws(rownames(pca_post))

# -----------------------------
# Load annotation
# -----------------------------
annot <- fread(ANNOT_CSV, data.table = FALSE)
colnames(annot) <- trimws(colnames(annot))

keep_ids <- readLines(KEEP_IDS_TXT)
keep_ids <- trimws(keep_ids)
keep_ids <- keep_ids[keep_ids != "Study.ID"]

annot$`Study ID` <- trimws(annot$`Study ID`)
annot$`Epic_Complete Bar code` <- trimws(annot$`Epic_Complete Bar code`)

annot_sub <- annot[annot$`Study ID` %in% keep_ids, , drop = FALSE]

barcodes <- annot_sub$`Epic_Complete Bar code`

covariates <- as.data.frame(annot_sub)
rownames(covariates) <- barcodes

# -----------------------------
# Align samples
# -----------------------------
shared <- Reduce(intersect, list(
  rownames(pca_pre),
  rownames(pca_post),
  rownames(covariates)
))

pca_pre    <- pca_pre[shared, , drop = FALSE]
pca_post   <- pca_post[shared, , drop = FALSE]
covariates <- covariates[shared, , drop = FALSE]

stopifnot(identical(rownames(pca_pre), rownames(covariates)))
stopifnot(identical(rownames(pca_post), rownames(covariates)))

# -----------------------------
# Check PC columns
# -----------------------------
message("PRE PCA dimensions:")
print(dim(pca_pre))

message("POST PCA dimensions:")
print(dim(pca_post))

if (ncol(pca_pre) < N_PCS || ncol(pca_post) < N_PCS) {
  stop("Not enough PC columns for N_PCS = ", N_PCS)
}

# -----------------------------
# Select PCs
# -----------------------------
pcs_pre  <- as.matrix(pca_pre[, 1:N_PCS, drop = FALSE])
pcs_post <- as.matrix(pca_post[, 1:N_PCS, drop = FALSE])

# -----------------------------
# Run UMAP
# -----------------------------
set.seed(42)

umap_pre <- umap(
  pcs_pre,
  n_neighbors = 15,
  min_dist = 0.1,
  metric = "euclidean",
  verbose = TRUE
)

umap_post <- umap(
  pcs_post,
  n_neighbors = 15,
  min_dist = 0.1,
  metric = "euclidean",
  verbose = TRUE
)

rownames(umap_pre) <- rownames(pcs_pre)
rownames(umap_post) <- rownames(pcs_post)

colnames(umap_pre) <- c("UMAP1", "UMAP2")
colnames(umap_post) <- c("UMAP1", "UMAP2")

# -----------------------------
# Collapse cycle phase: PE vs Secretory only
# -----------------------------
phase <- trimws(as.character(covariates$`Cycle phase for Analysis`))

covariates$Cycle_binary <- NA_character_
covariates$Cycle_binary[phase == "PE"] <- "Proliferative"
covariates$Cycle_binary[phase %in% c("ESE", "MSE", "LSE", "SE")] <- "Secretory"

keep_cycle <- !is.na(covariates$Cycle_binary)

covariates_cycle <- covariates[keep_cycle, , drop = FALSE]
umap_pre_cycle   <- umap_pre[keep_cycle, , drop = FALSE]
umap_post_cycle  <- umap_post[keep_cycle, , drop = FALSE]

covariates_cycle$Cycle_binary <- factor(
  covariates_cycle$Cycle_binary,
  levels = c("Proliferative", "Secretory")
)

# -----------------------------
# Plot helper: general UMAP
# -----------------------------
plot_umap <- function(coords, labels, title, file) {
  
  labs <- as.factor(labels)
  cols <- rainbow(length(levels(labs)))[labs]
  
  png(file, width = 2200, height = 1600, res = 300)
  
  plot(
    coords[, 1],
    coords[, 2],
    pch = 16,
    col = cols,
    xlab = "UMAP1",
    ylab = "UMAP2",
    main = title
  )
  
  legend(
    "topright",
    legend = levels(labs),
    col = rainbow(length(levels(labs))),
    pch = 16,
    cex = 0.8
  )
  
  dev.off()
}

# -----------------------------
# Plot helper: case/control UMAP
# -----------------------------
plot_umap_case_control <- function(coords, labels, title, file) {
  
  labs <- as.character(labels)
  cols <- ifelse(labs == "Yes", "red", "black")
  
  png(file, width = 2200, height = 1600, res = 300)
  
  plot(
    coords[, 1],
    coords[, 2],
    pch = 16,
    col = cols,
    xlab = "UMAP1",
    ylab = "UMAP2",
    main = title
  )
  
  legend(
    "topright",
    legend = c("No", "Yes"),
    col = c("black", "red"),
    pch = 16,
    cex = 0.9
  )
  
  dev.off()
}

# -----------------------------
# Plot helper: binary cycle phase UMAP
# -----------------------------
plot_umap_binary_cycle <- function(coords, labels, title, file) {
  
  labs <- factor(labels, levels = c("Proliferative", "Secretory"))
  
  cycle_colors <- c(
    "Proliferative" = "blue",
    "Secretory" = "orange"
  )
  
  cols <- cycle_colors[as.character(labs)]
  
  png(file, width = 2200, height = 1600, res = 300)
  
  plot(
    coords[, 1],
    coords[, 2],
    pch = 16,
    col = cols,
    xlab = "UMAP1",
    ylab = "UMAP2",
    main = title
  )
  
  legend(
    "topright",
    legend = levels(labs),
    col = cycle_colors[levels(labs)],
    pch = 16,
    cex = 0.9
  )
  
  dev.off()
}

# -----------------------------
# Save plots: case/control
# -----------------------------
plot_umap_case_control(
  umap_pre,
  covariates$`Endometriosis (Yes/No)`,
  "UMAP (PRE SmartSVA): Case vs Control",
  file.path(OUT_DIR, "UMAP_pre_case_control.png")
)

plot_umap_case_control(
  umap_post,
  covariates$`Endometriosis (Yes/No)`,
  "UMAP (POST SmartSVA): Case vs Control",
  file.path(OUT_DIR, "UMAP_post_case_control.png")
)

# -----------------------------
# Save plots: original multi-level cycle phase
# -----------------------------
plot_umap(
  umap_pre,
  covariates$`Cycle phase for Analysis`,
  "UMAP (PRE SmartSVA): Cycle Phase",
  file.path(OUT_DIR, "UMAP_pre_cycle_phase.png")
)

plot_umap(
  umap_post,
  covariates$`Cycle phase for Analysis`,
  "UMAP (POST SmartSVA): Cycle Phase",
  file.path(OUT_DIR, "UMAP_post_cycle_phase.png")
)

# -----------------------------
# Save plots: PE vs Secretory only
# -----------------------------
plot_umap_binary_cycle(
  umap_pre_cycle,
  covariates_cycle$Cycle_binary,
  "UMAP (PRE SmartSVA): Proliferative vs Secretory",
  file.path(OUT_DIR_BINARY, "UMAP_pre_proliferative_vs_secretory_no_menstrual.png")
)

plot_umap_binary_cycle(
  umap_post_cycle,
  covariates_cycle$Cycle_binary,
  "UMAP (POST SmartSVA): Proliferative vs Secretory",
  file.path(OUT_DIR_BINARY, "UMAP_post_proliferative_vs_secretory_no_menstrual.png")
)

# -----------------------------
# Save plots: batch
# -----------------------------
plot_umap(
  umap_pre,
  covariates$`Batch`,
  "UMAP (PRE SmartSVA): Batch",
  file.path(OUT_DIR, "UMAP_pre_batch.png")
)

plot_umap(
  umap_post,
  covariates$`Batch`,
  "UMAP (POST SmartSVA): Batch",
  file.path(OUT_DIR, "UMAP_post_batch.png")
)

# -----------------------------
# Save plots: institute for analysis
# -----------------------------
plot_umap(
  umap_pre,
  covariates$`Institute for Analysis`,
  "UMAP (PRE SmartSVA): Institute for Analysis",
  file.path(OUT_DIR, "UMAP_pre_institute_for_analysis.png")
)

plot_umap(
  umap_post,
  covariates$`Institute for Analysis`,
  "UMAP (POST SmartSVA): Institute for Analysis",
  file.path(OUT_DIR, "UMAP_post_institute_for_analysis.png")
)

# -----------------------------
# Save UMAP coordinates
# -----------------------------
umap_pre_df <- as.data.frame(umap_pre)
umap_pre_df$Sample <- rownames(umap_pre)

write.csv(
  umap_pre_df,
  file.path(OUT_DIR, "umap_pre_coordinates.csv"),
  row.names = FALSE
)

umap_post_df <- as.data.frame(umap_post)
umap_post_df$Sample <- rownames(umap_post)

write.csv(
  umap_post_df,
  file.path(OUT_DIR, "umap_post_coordinates.csv"),
  row.names = FALSE
)

# -----------------------------
# Save PE vs Secretory-only UMAP coordinates
# -----------------------------
umap_pre_cycle_df <- as.data.frame(umap_pre_cycle)
umap_pre_cycle_df$Sample <- rownames(umap_pre_cycle)
umap_pre_cycle_df$Cycle_binary <- covariates_cycle$Cycle_binary

write.csv(
  umap_pre_cycle_df,
  file.path(OUT_DIR_BINARY, "umap_pre_PE_SE_only_coordinates.csv"),
  row.names = FALSE
)

umap_post_cycle_df <- as.data.frame(umap_post_cycle)
umap_post_cycle_df$Sample <- rownames(umap_post_cycle)
umap_post_cycle_df$Cycle_binary <- covariates_cycle$Cycle_binary

write.csv(
  umap_post_cycle_df,
  file.path(OUT_DIR_BINARY, "umap_post_PE_SE_only_coordinates.csv"),
  row.names = FALSE
)

# -----------------------------
# Final checks
# -----------------------------
cat("\nOriginal cycle phase counts:\n")
print(table(phase, useNA = "ifany"))

cat("\nIncluded Proliferative vs Secretory counts:\n")
print(table(covariates_cycle$Cycle_binary, useNA = "ifany"))

cat("\nExcluded phase counts:\n")
print(table(phase[!keep_cycle], useNA = "ifany"))

cat("\nMain UMAP files written to:\n")
print(OUT_DIR)
print(list.files(OUT_DIR))

cat("\nPE vs Secretory-only UMAP files written to:\n")
print(OUT_DIR_BINARY)
print(list.files(OUT_DIR_BINARY))