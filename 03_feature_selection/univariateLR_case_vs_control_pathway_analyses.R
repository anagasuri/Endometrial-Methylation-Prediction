# -----------------------------
# CASE VS CONTROL PATHWAY ANALYSIS
# Univariate LR FDR-significant CpGs only
# GO / KEGG / Reactome via missMethyl
# -----------------------------

library(missMethyl)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(reactome.db)
library(dplyr)
library(ggplot2)

# -----------------------------
# Input files
# -----------------------------

sig_file <- "/mnt/efs/home/ubuntu/endo-me_data/ttest_feature_selection_outputs/case_control_univariate_logistic/case_control_univariate_logistic_FDR_0.05_selected_cpgs_n2.csv"

all_file <- "/mnt/efs/home/ubuntu/endo-me_data/ttest_feature_selection_outputs/case_control_univariate_logistic/case_control_univariate_logistic_all_cpgs_with_FDR.csv"

out_dir <- "/mnt/efs/home/ubuntu/endo-me_data/ttest_feature_selection_outputs/ttest_pathway_analysis_outputs/case_control/univariate_LR_FDR"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n====================================================\n")
cat("Starting case/control univariate LR pathway analysis\n")
cat("Start time:", as.character(Sys.time()), "\n")
cat("Output directory:", out_dir, "\n")
cat("====================================================\n\n")

# -----------------------------
# Helper: get CpG column
# -----------------------------

get_cpg_vector <- function(df) {
  possible_cols <- c("cpg", "CpG", "probe", "Probe", "CpG_ID", "cpg_id")
  found <- possible_cols[possible_cols %in% colnames(df)]
  
  if (length(found) == 0) {
    stop(
      paste(
        "Could not find CpG column. Found columns:",
        paste(colnames(df), collapse = ", ")
      )
    )
  }
  
  cat("Using CpG column:", found[1], "\n")
  return(unique(na.omit(df[[found[1]]])))
}

# -----------------------------
# Helper: save FDR significant pathways
# -----------------------------

save_fdr_results <- function(df, db_name, out_file) {
  fdr_df <- df %>%
    filter(FDR < 0.05) %>%
    arrange(FDR, P.DE)
  
  write.csv(fdr_df, out_file, row.names = FALSE)
  
  cat("\n", db_name, "FDR-significant pathways:", nrow(fdr_df), "\n")
  cat("Saved:", out_file, "\n")
  
  if (nrow(fdr_df) > 0) {
    print(fdr_df)
  } else {
    cat("No FDR-significant", db_name, "pathways found.\n")
  }
  
  return(fdr_df)
}

# -----------------------------
# Helper: plot FDR significant pathways only
# -----------------------------

plot_fdr_pathways <- function(df, db_name, out_file) {
  
  plot_df <- df %>%
    filter(FDR < 0.05) %>%
    arrange(FDR) %>%
    head(20)
  
  if (nrow(plot_df) == 0) {
    cat("Skipping", db_name, "plot because no pathways passed FDR < 0.05.\n")
    return(NULL)
  }
  
  if (db_name == "GO") {
    y_col <- "TERM"
  } else if (db_name == "KEGG") {
    y_col <- "Description"
  } else if (db_name == "Reactome") {
    if (!"Pathway" %in% colnames(plot_df)) {
      plot_df <- plot_df %>% mutate(Pathway = rownames(plot_df))
    }
    y_col <- "Pathway"
  }
  
  p <- ggplot(
    plot_df,
    aes(
      x = -log10(FDR),
      y = reorder(.data[[y_col]], -log10(FDR))
    )
  ) +
    geom_bar(stat = "identity") +
    labs(
      title = paste0(db_name, " FDR-significant pathways"),
      x = "-log10(FDR)",
      y = ""
    ) +
    theme_minimal()
  
  print(p)
  ggsave(out_file, plot = p, width = 10, height = 7, dpi = 300)
  
  cat("Saved plot:", out_file, "\n")
  
  return(p)
}

# -----------------------------
# Load CpG files
# -----------------------------

cat("Loading CpG files...\n")
cat("Time:", as.character(Sys.time()), "\n\n")

sig_df <- read.csv(sig_file, stringsAsFactors = FALSE)
all_df <- read.csv(all_file, stringsAsFactors = FALSE)

cat("Selected CpG file dimensions:", dim(sig_df), "\n")
cat("Background CpG file dimensions:", dim(all_df), "\n\n")

sig.cpg <- get_cpg_vector(sig_df)
all.cpg <- get_cpg_vector(all_df)

cat("\nSelected FDR CpGs:", length(sig.cpg), "\n")
cat("Background CpGs:", length(all.cpg), "\n")
cat("All selected CpGs in background:", all(sig.cpg %in% all.cpg), "\n")

missing_cpgs <- setdiff(sig.cpg, all.cpg)

if (length(missing_cpgs) > 0) {
  cat("Warning:", length(missing_cpgs), "selected CpGs are missing from background.\n")
  cat("Keeping only selected CpGs that exist in background.\n")
  sig.cpg <- intersect(sig.cpg, all.cpg)
}

cat("Final selected CpGs used:", length(sig.cpg), "\n\n")

# -----------------------------
# CpG to Entrez mapping
# -----------------------------

cat("Mapping CpGs to Entrez IDs...\n")
cat("Start time:", as.character(Sys.time()), "\n")

mapped <- getMappedEntrezIDs(
  sig.cpg = sig.cpg,
  all.cpg = all.cpg,
  array.type = "EPIC"
)

sig_entrez <- unique(mapped$sig.eg)
sig_entrez <- sig_entrez[!is.na(sig_entrez)]

cat("Finished CpG to Entrez mapping.\n")
cat("End time:", as.character(Sys.time()), "\n")
cat("Mapped significant Entrez genes:", length(sig_entrez), "\n\n")

gene_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = as.character(sig_entrez),
  keytype = "ENTREZID",
  columns = c("ENTREZID", "SYMBOL", "GENENAME")
) %>%
  distinct()

write.csv(
  gene_map,
  file.path(out_dir, "case_control_FDR_univariate_LR_mapped_entrez_genes.csv"),
  row.names = FALSE
)

cat("Mapped genes saved.\n\n")

# -----------------------------
# GO enrichment
# -----------------------------

cat("Running GO enrichment...\n")
cat("Start time:", as.character(Sys.time()), "\n")

go_res <- gometh(
  sig.cpg = sig.cpg,
  all.cpg = all.cpg,
  collection = "GO",
  array.type = "EPIC",
  prior.prob = TRUE,
  fract.counts = TRUE,
  plot.bias = FALSE,
  sig.genes = TRUE
)

cat("Finished GO enrichment.\n")
cat("End time:", as.character(Sys.time()), "\n\n")

write.csv(
  go_res,
  file.path(out_dir, "case_control_FDR_univariate_LR_GO_all_results.csv"),
  row.names = FALSE
)

cat("Saved GO all results.\n\n")

# -----------------------------
# KEGG enrichment
# -----------------------------

cat("Running KEGG enrichment...\n")
cat("Start time:", as.character(Sys.time()), "\n")

kegg_res <- gometh(
  sig.cpg = sig.cpg,
  all.cpg = all.cpg,
  collection = "KEGG",
  array.type = "EPIC",
  prior.prob = TRUE,
  fract.counts = TRUE,
  plot.bias = FALSE,
  sig.genes = TRUE
)

cat("Finished KEGG enrichment.\n")
cat("End time:", as.character(Sys.time()), "\n\n")

write.csv(
  kegg_res,
  file.path(out_dir, "case_control_FDR_univariate_LR_KEGG_all_results.csv"),
  row.names = FALSE
)

cat("Saved KEGG all results.\n\n")

# -----------------------------
# Build Reactome gene sets
# -----------------------------

cat("Building Reactome gene sets...\n")
cat("Start time:", as.character(Sys.time()), "\n")

reactome_sets <- AnnotationDbi::as.list(reactomeEXTID2PATHID)

reactome_pathway2gene <- split(
  rep(names(reactome_sets), lengths(reactome_sets)),
  unlist(reactome_sets)
)

reactome_names <- AnnotationDbi::as.list(reactomePATHID2NAME)

reactome_pathway2gene <- reactome_pathway2gene[
  names(reactome_pathway2gene) %in% names(reactome_names)
]

names(reactome_pathway2gene) <- reactome_names[names(reactome_pathway2gene)]

reactome_pathway2gene <- reactome_pathway2gene[
  lengths(reactome_pathway2gene) >= 10 & lengths(reactome_pathway2gene) <= 500
]

cat("Finished building Reactome gene sets.\n")
cat("End time:", as.character(Sys.time()), "\n")
cat("Reactome pathways retained:", length(reactome_pathway2gene), "\n\n")

# -----------------------------
# Reactome enrichment
# -----------------------------

cat("Running Reactome enrichment...\n")
cat("Start time:", as.character(Sys.time()), "\n")

reactome_res <- gsameth(
  sig.cpg = sig.cpg,
  all.cpg = all.cpg,
  collection = reactome_pathway2gene,
  array.type = "EPIC",
  prior.prob = TRUE,
  fract.counts = TRUE,
  plot.bias = FALSE,
  sig.genes = TRUE
)

cat("Finished Reactome enrichment.\n")
cat("End time:", as.character(Sys.time()), "\n\n")

reactome_res <- reactome_res %>%
  mutate(Pathway = rownames(.))

write.csv(
  reactome_res,
  file.path(out_dir, "case_control_FDR_univariate_LR_Reactome_all_results.csv"),
  row.names = FALSE
)

cat("Saved Reactome all results.\n\n")

# -----------------------------
# Save FDR-significant pathway tables only
# -----------------------------

cat("Saving FDR-significant pathway tables...\n")

go_fdr <- save_fdr_results(
  go_res,
  db_name = "GO",
  out_file = file.path(out_dir, "case_control_FDR_univariate_LR_GO_FDR_significant.csv")
)

kegg_fdr <- save_fdr_results(
  kegg_res,
  db_name = "KEGG",
  out_file = file.path(out_dir, "case_control_FDR_univariate_LR_KEGG_FDR_significant.csv")
)

reactome_fdr <- save_fdr_results(
  reactome_res,
  db_name = "Reactome",
  out_file = file.path(out_dir, "case_control_FDR_univariate_LR_Reactome_FDR_significant.csv")
)

# -----------------------------
# Plot FDR-significant pathways only
# -----------------------------

cat("\nPlotting FDR-significant pathways if any exist...\n")

plot_fdr_pathways(
  go_res,
  db_name = "GO",
  out_file = file.path(out_dir, "case_control_FDR_univariate_LR_GO_FDR_barplot.png")
)

plot_fdr_pathways(
  kegg_res,
  db_name = "KEGG",
  out_file = file.path(out_dir, "case_control_FDR_univariate_LR_KEGG_FDR_barplot.png")
)

plot_fdr_pathways(
  reactome_res,
  db_name = "Reactome",
  out_file = file.path(out_dir, "case_control_FDR_univariate_LR_Reactome_FDR_barplot.png")
)

# -----------------------------
# Final summary
# -----------------------------

cat("\n====================================================\n")
cat("Case/control univariate LR FDR CpG pathway analysis complete.\n")
cat("End time:", as.character(Sys.time()), "\n")
cat("Output directory:\n")
cat(out_dir, "\n\n")

cat("Selected FDR CpGs:", length(sig.cpg), "\n")
cat("Mapped significant Entrez genes:", length(sig_entrez), "\n")
cat("GO FDR-significant pathways:", nrow(go_fdr), "\n")
cat("KEGG FDR-significant pathways:", nrow(kegg_fdr), "\n")
cat("Reactome FDR-significant pathways:", nrow(reactome_fdr), "\n")
cat("====================================================\n")