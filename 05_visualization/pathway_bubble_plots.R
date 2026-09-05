############################################################
# PATHWAY BUBBLE PLOTS
# Cycle phase = FDR-significant pathways
# Case/control = nominal pathways only
# Saves to: /mnt/efs/home/ubuntu/defense_final_figs
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
})

out_dir <- "/mnt/efs/home/ubuntu/defense_final_figs"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

############################################################
# FILE PATHS
############################################################

files <- list(
  case_control_ttest_GO = "/mnt/efs/home/ubuntu/endo-me_data/ttest_feature_selection_outputs/ttest_pathway_analysis_outputs/case_control/FDR_ttest/case_control_FDR_ttest_GO_all_results.csv",
  case_control_ttest_KEGG = "/mnt/efs/home/ubuntu/endo-me_data/ttest_feature_selection_outputs/ttest_pathway_analysis_outputs/case_control/FDR_ttest/case_control_FDR_ttest_KEGG_all_results.csv",
  case_control_ttest_Reactome = "/mnt/efs/home/ubuntu/endo-me_data/ttest_feature_selection_outputs/ttest_pathway_analysis_outputs/case_control/FDR_ttest/case_control_FDR_ttest_Reactome_all_results.csv",
  
  cycle_phase_ttest_GO = "/mnt/efs/home/ubuntu/endo-me_data/ttest_feature_selection_outputs/ttest_pathway_analysis_outputs/cycle_phase/FDR_ttest/cycle_phase_FDR_ttest_GO_all_results.csv",
  cycle_phase_ttest_KEGG = "/mnt/efs/home/ubuntu/endo-me_data/ttest_feature_selection_outputs/ttest_pathway_analysis_outputs/cycle_phase/FDR_ttest/cycle_phase_FDR_ttest_KEGG_all_results.csv",
  cycle_phase_ttest_Reactome = "/mnt/efs/home/ubuntu/endo-me_data/ttest_feature_selection_outputs/ttest_pathway_analysis_outputs/cycle_phase/FDR_ttest/cycle_phase_FDR_ttest_Reactome_all_results.csv",
  
  case_control_univLR_GO = "/mnt/efs/home/ubuntu/endo-me_data/ttest_feature_selection_outputs/ttest_pathway_analysis_outputs/case_control/univariate_LR_FDR/case_control_FDR_univariate_LR_GO_all_results.csv",
  case_control_univLR_KEGG = "/mnt/efs/home/ubuntu/endo-me_data/ttest_feature_selection_outputs/ttest_pathway_analysis_outputs/case_control/univariate_LR_FDR/case_control_FDR_univariate_LR_KEGG_all_results.csv",
  case_control_univLR_Reactome = "/mnt/efs/home/ubuntu/endo-me_data/ttest_feature_selection_outputs/ttest_pathway_analysis_outputs/case_control/univariate_LR_FDR/case_control_FDR_univariate_LR_Reactome_all_results.csv",
  
  # NOTE: mislabeled earlier — using ttest files for univLR
  cycle_phase_univLR_GO = "/mnt/efs/home/ubuntu/endo-me_data/ttest_feature_selection_outputs/ttest_pathway_analysis_outputs/cycle_phase/univariate_LR_FDR/cycle_phase_FDR_ttest_GO_all_results.csv",
  cycle_phase_univLR_KEGG = "/mnt/efs/home/ubuntu/endo-me_data/ttest_feature_selection_outputs/ttest_pathway_analysis_outputs/cycle_phase/univariate_LR_FDR/cycle_phase_FDR_ttest_KEGG_all_results.csv",
  cycle_phase_univLR_Reactome = "/mnt/efs/home/ubuntu/endo-me_data/ttest_feature_selection_outputs/ttest_pathway_analysis_outputs/cycle_phase/univariate_LR_FDR/cycle_phase_FDR_ttest_Reactome_all_results.csv"
)

############################################################
# SHARED FUNCTION
############################################################

make_bubble_plot <- function(file, label, database, mode = c("FDR", "Nominal"), top_n = 20) {
  
  mode <- match.arg(mode)
  
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
  
  if (mode == "FDR") {
    plot_df <- df %>%
      filter(!is.na(FDR), FDR < 0.05) %>%
      arrange(FDR, P.DE) %>%
      slice_head(n = top_n) %>%
      mutate(
        pathway_name = str_wrap(.data[[pathway_col]], 45),
        x_value = -log10(FDR),
        color_value = FDR,
        gene_count = .data[[size_col]]
      )
    
    x_lab <- "-log10(FDR)"
    color_lab <- "FDR"
    title_text <- paste(label, database, "FDR-Significant Pathway Enrichment")
    suffix <- "FDR_bubble"
    
  } else {
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
    
    x_lab <- "-log10(P.DE)"
    color_lab <- "P.DE"
    title_text <- paste(label, database, "Nominal Pathway Enrichment")
    suffix <- "nominal_bubble"
  }
  
  if (nrow(plot_df) == 0) {
    message("No pathways for: ", label, " ", database, " using mode: ", mode)
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
      title = title_text,
      x = x_lab,
      y = "",
      size = "Genes",
      color = color_lab
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.y = element_text(size = 9)
    )
  
  out_file <- file.path(out_dir, paste0(label, "_", database, "_", suffix, ".png"))
  ggsave(out_file, plot = p, width = 10, height = 7, dpi = 300)
  
  message("Saved: ", out_file)
  return(p)
}

############################################################
# RUN CASE/CONTROL AS NOMINAL
############################################################

make_bubble_plot(files$case_control_ttest_GO, "case_control_ttest", "GO", mode = "Nominal")
make_bubble_plot(files$case_control_ttest_KEGG, "case_control_ttest", "KEGG", mode = "Nominal")
make_bubble_plot(files$case_control_ttest_Reactome, "case_control_ttest", "Reactome", mode = "Nominal")

make_bubble_plot(files$case_control_univLR_GO, "case_control_univLR", "GO", mode = "Nominal")
make_bubble_plot(files$case_control_univLR_KEGG, "case_control_univLR", "KEGG", mode = "Nominal")
make_bubble_plot(files$case_control_univLR_Reactome, "case_control_univLR", "Reactome", mode = "Nominal")

############################################################
# RUN CYCLE PHASE AS FDR-SIGNIFICANT
############################################################

make_bubble_plot(files$cycle_phase_ttest_GO, "cycle_phase_ttest", "GO", mode = "FDR")
make_bubble_plot(files$cycle_phase_ttest_KEGG, "cycle_phase_ttest", "KEGG", mode = "FDR")
make_bubble_plot(files$cycle_phase_ttest_Reactome, "cycle_phase_ttest", "Reactome", mode = "FDR")

make_bubble_plot(files$cycle_phase_univLR_GO, "cycle_phase_univLR", "GO", mode = "FDR")
make_bubble_plot(files$cycle_phase_univLR_KEGG, "cycle_phase_univLR", "KEGG", mode = "FDR")
make_bubble_plot(files$cycle_phase_univLR_Reactome, "cycle_phase_univLR", "Reactome", mode = "FDR")

cat("\nDONE. Bubble plots saved to:\n", out_dir, "\n")