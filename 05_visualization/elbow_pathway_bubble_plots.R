############################################################
# ELBOW PATHWAY BUBBLE PLOTS
# Case/control + Cycle phase
# All plots = nominal pathways only (P.DE < 0.05)
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
})

############################################################
# OUTPUT DIRECTORY
############################################################

out_dir <- "/mnt/efs/home/ubuntu/defense_final_figs/elbow_nominal_bubbleplots"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

############################################################
# FILE PATHS
############################################################

files <- list(
  case_control_elbow_GO =
    "/mnt/efs/home/ubuntu/case_control_top978_GO_analysis/case_control_top978_GO_enrichment_all.csv",
  
  case_control_elbow_KEGG =
    "/mnt/efs/home/ubuntu/case_control_top978_KEGG_analysis/case_control_top978_KEGG_enrichment_all.csv",
  
  cycle_phase_elbow_GO =
    "/mnt/efs/home/ubuntu/cycle_phase_top2303_GO_analysis/cycle_phase_top2303_GO_enrichment_all.csv",
  
  cycle_phase_elbow_KEGG =
    "/mnt/efs/home/ubuntu/cycle_phase_top2303_KEGG_analysis/cycle_phase_top2303_KEGG_enrichment_all.csv"
)

############################################################
# SHARED FUNCTION
############################################################

make_bubble_plot <- function(file, label, database, top_n = 20) {
  
  if (!file.exists(file)) {
    message("File does not exist: ", file)
    return(NULL)
  }
  
  df <- read_csv(file, show_col_types = FALSE)
  
  pathway_col <- case_when(
    "TERM" %in% colnames(df) ~ "TERM",
    "Description" %in% colnames(df) ~ "Description",
    "Pathway" %in% colnames(df) ~ "Pathway",
    TRUE ~ NA_character_
  )
  
  size_col <- case_when(
    "DE" %in% colnames(df) ~ "DE",
    "SigGenesInSet" %in% colnames(df) ~ "SigGenesInSet",
    "N" %in% colnames(df) ~ "N",
    TRUE ~ NA_character_
  )
  
  if (is.na(pathway_col) | is.na(size_col)) {
    message("Missing pathway or size column in: ", file)
    return(NULL)
  }
  
  plot_df <- df %>%
    filter(!is.na(P.DE), P.DE < 0.05) %>%
    arrange(P.DE) %>%
    slice_head(n = top_n) %>%
    mutate(
      pathway_name = str_wrap(.data[[pathway_col]], 45),
      x_value = -log10(P.DE),
      color_value = P.DE,
      gene_count = .data[[size_col]]
    )
  
  if (nrow(plot_df) == 0) {
    message("No nominal pathways for: ", label, " ", database)
    return(NULL)
  }
  
  p <- ggplot(
    plot_df,
    aes(
      x = x_value,
      y = reorder(pathway_name, x_value),
      size = gene_count,
      color = color_value
    )
  ) +
    geom_point(alpha = 0.85) +
    scale_color_gradient(low = "red", high = "blue") +
    labs(
      title = paste(label, database, "Nominal Pathway Enrichment"),
      x = "-log10(P.DE)",
      y = "",
      size = "Genes",
      color = "P.DE"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.y = element_text(size = 9)
    )
  
  out_file <- file.path(
    out_dir,
    paste0(label, "_", database, "_nominal_bubbleplot.png")
  )
  
  ggsave(out_file, plot = p, width = 10, height = 7, dpi = 300)
  
  message("Saved: ", out_file)
  return(p)
}

############################################################
# RUN CASE/CONTROL ELBOW NOMINAL PLOTS
############################################################

make_bubble_plot(
  files$case_control_elbow_GO,
  "Case_Control_Elbow_Top978",
  "GO",
  top_n = 20
)

make_bubble_plot(
  files$case_control_elbow_KEGG,
  "Case_Control_Elbow_Top978",
  "KEGG",
  top_n = 20
)

############################################################
# RUN CYCLE PHASE ELBOW NOMINAL PLOTS
############################################################

make_bubble_plot(
  files$cycle_phase_elbow_GO,
  "Cycle_Phase_Elbow_Top2303",
  "GO",
  top_n = 20
)

make_bubble_plot(
  files$cycle_phase_elbow_KEGG,
  "Cycle_Phase_Elbow_Top2303",
  "KEGG",
  top_n = 20
)

cat("\nDONE. Nominal elbow bubble plots saved to:\n", out_dir, "\n")