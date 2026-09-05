############################################################
# Heatmap of Top 978 Predictive CpGs - Case/Control (ELBOW)
############################################################

suppressPackageStartupMessages({
  library(reticulate)
  library(pheatmap)
  library(dplyr)
})

############################################################
# FILE PATHS 
############################################################

RANKED_CPG_FILE <- "/mnt/efs/home/ubuntu/case_control_top_cpgs/case_control_ranked_cpg_coefficients.csv"
H5AD_FILE <- "/mnt/efs/home/ubuntu/M_smartsva_corrected_with_SH_annot_filtered.h5ad"
OUT_FILE <- "/mnt/efs/home/ubuntu/defense_final_figs/case_control_top978_heatmap_unforced.png"

N_CPGS <- 978

############################################################
# Load Python modules
############################################################

ad <- import("anndata", convert = FALSE)
scipy_sparse <- import("scipy.sparse", convert = FALSE)

############################################################
# Read ranked CpG file
############################################################

ranked_df <- read.csv(RANKED_CPG_FILE, stringsAsFactors = FALSE)

candidate_cols <- c("cpg", "CpG", "cpg_id", "probe", "probe_id", "ilmn", "Name")
cpg_col <- intersect(candidate_cols, colnames(ranked_df))[1]

if (is.na(cpg_col)) stop("No CpG column found.")

if ("rank" %in% colnames(ranked_df)) {
  ranked_df <- ranked_df[order(ranked_df$rank), ]
} else {
  ranked_df <- ranked_df[order(ranked_df$rank_by_abs_coefficient), ]
}

selected_cpgs <- unique(na.omit(ranked_df[[cpg_col]]))[1:N_CPGS]

############################################################
# Read h5ad
############################################################

adata <- ad$read_h5ad(H5AD_FILE)

meta <- py_to_r(adata$obs)
meta$sample_id <- rownames(meta)

cpg_names <- py_to_r(adata$var_names$to_list())
sample_names <- py_to_r(adata$obs_names$to_list())

X <- adata$X

if (py_to_r(scipy_sparse$issparse(X))) {
  meth <- py_to_r(X$toarray())
} else {
  meth <- py_to_r(X)
}

meth <- as.matrix(meth)
rownames(meth) <- sample_names
colnames(meth) <- cpg_names

# transpose → CpGs x samples
meth <- t(meth)

############################################################
# Subset to selected CpGs present in matrix
############################################################

selected_cpgs_present <- selected_cpgs[selected_cpgs %in% rownames(meth)]

if (length(selected_cpgs_present) == 0) {
  stop("None of the selected CpGs were found in the h5ad matrix.")
}

meth_sub <- meth[selected_cpgs_present, , drop = FALSE]

cat("Number of CpGs plotted:", nrow(meth_sub), "\n")
cat("Number of samples plotted:", ncol(meth_sub), "\n")

############################################################
# Metadata labels
############################################################

case_col <- intersect(c("endo", "case_control", "group"), colnames(meta))[1]

if (is.na(case_col)) stop("No case/control column found.")

phase_col <- intersect(
  c("cycle-phase", "CyclePhase", "phase", "Phase", "cycle", "Cycle"),
  colnames(meta)
)[1]

if (is.na(phase_col)) stop("No cycle phase column found.")

meta <- meta[match(colnames(meth_sub), meta$sample_id), , drop = FALSE]

if (any(is.na(meta$sample_id))) {
  stop("Some samples in matrix were not found in metadata.")
}

meta$CaseControl <- factor(meta[[case_col]], levels = c(0, 1))

meta$CyclePhase <- dplyr::case_when(
  meta[[phase_col]] %in% c(
    "PE", "P", "Proliferative", "proliferative",
    "Early proliferative", "Late proliferative"
  ) ~ "Proliferative",
  meta[[phase_col]] %in% c(
    "ESE", "MSE", "LSE", "SE", "Secretory", "secretory",
    "Early secretory", "Mid secretory", "Late secretory"
  ) ~ "Secretory",
  TRUE ~ NA_character_
)

annotation_col <- data.frame(
  CyclePhase = factor(meta$CyclePhase, levels = c("Proliferative", "Secretory")),
  CaseControl = meta$CaseControl
)

rownames(annotation_col) <- meta$sample_id

############################################################
# Z-score CpGs across samples
############################################################

mode(meth_sub) <- "numeric"

meth_scaled <- t(scale(t(meth_sub)))
meth_scaled[is.na(meth_scaled)] <- 0

############################################################
# Colors
############################################################

heat_colors <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
breaks <- seq(-3, 3, length.out = 101)

ann_colors <- list(
  CyclePhase = c(
    "Proliferative" = "blue",
    "Secretory" = "magenta"
  ),
  CaseControl = c(
    "0" = "black",
    "1" = "red"
  )
)

############################################################
# Save heatmap
############################################################

pheatmap(
  meth_scaled,
  color = heat_colors,
  breaks = breaks,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  show_rownames = FALSE,
  show_colnames = FALSE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  treeheight_row = 0,
  fontsize = 12,
  main = "Top 978 Predictive CpGs - Case vs Control",
  filename = OUT_FILE,
  width = 20,
  height = 18
)

cat("Saved to:", OUT_FILE, "\n")
cat("File exists:", file.exists(OUT_FILE), "\n")