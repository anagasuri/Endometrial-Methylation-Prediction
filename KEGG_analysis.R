############################################################
# CpG -> Entrez gene mapping + KEGG pathway analysis
# for EPIC methylation array data
############################################################

# If needed:
# install.packages("BiocManager")
# BiocManager::install("missMethyl", ask = FALSE)
# BiocManager::install("org.Hs.eg.db", ask = FALSE)

suppressPackageStartupMessages({
  library(missMethyl)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

############################################################
# Input files
############################################################

TOP_CPG_CSV <- "/home/ubuntu/ridge_top_1000_features.csv"
ALL_CPG_CSV <- "/home/ubuntu/ridge_top_features_all.csv"

OUT_PREFIX <- "ridge_top1000"
OUT_DIR <- "KEGG_analysis"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

############################################################
# Read input files
############################################################

top_df <- read.csv(TOP_CPG_CSV, stringsAsFactors = FALSE)
all_df <- read.csv(ALL_CPG_CSV, stringsAsFactors = FALSE)

cat("Top CpG file dimensions:\n")
print(dim(top_df))
cat("All-feature file dimensions:\n")
print(dim(all_df))

cat("\nTop CpG file columns:\n")
print(colnames(top_df))

cat("\nAll-feature file columns:\n")
print(colnames(all_df))

############################################################
# Detect CpG column names
############################################################

candidate_cols <- c("cpg", "CpG", "cpg_id", "probe", "probe_id", "ilmn", "Name")

top_cpg_col <- intersect(candidate_cols, colnames(top_df))[1]
all_cpg_col <- intersect(candidate_cols, colnames(all_df))[1]

if (is.na(top_cpg_col) || is.na(all_cpg_col)) {
  stop("Could not find CpG column in one or both input files.")
}

cat("\nUsing CpG column in top file:\n")
print(top_cpg_col)

cat("\nUsing CpG column in all-feature file:\n")
print(all_cpg_col)

############################################################
# Extract CpG vectors
############################################################

top_cpgs <- unique(na.omit(top_df[[top_cpg_col]]))
all_cpgs <- unique(na.omit(all_df[[all_cpg_col]]))

cat("\nNumber of top CpGs:\n")
print(length(top_cpgs))

cat("\nNumber of all CpGs in universe:\n")
print(length(all_cpgs))

cat("\nFirst few top CpGs:\n")
print(head(top_cpgs))

############################################################
# Diagnostic checks
############################################################

cat("\nDiagnostic checks:\n")
cat("length(top_cpgs): ", length(top_cpgs), "\n")
cat("length(all_cpgs): ", length(all_cpgs), "\n")
cat("setequal(top_cpgs, all_cpgs): ", setequal(top_cpgs, all_cpgs), "\n")
cat("Top CpGs missing from universe: ", sum(!top_cpgs %in% all_cpgs), "\n")

if (setequal(top_cpgs, all_cpgs)) {
  warning("top_cpgs and all_cpgs are identical. This is usually not correct for enrichment testing.")
}

if (sum(!top_cpgs %in% all_cpgs) > 0) {
  warning("Some top CpGs are not present in all_cpgs. Check your universe file.")
}

############################################################
# Map top CpGs to Entrez IDs
############################################################

mapped <- getMappedEntrezIDs(
  sig.cpg = top_cpgs,
  all.cpg = all_cpgs,
  array.type = "EPIC"
)

entrez_ids <- unique(mapped$sig.eg)
entrez_ids <- entrez_ids[!is.na(entrez_ids)]

cat("\nNumber of mapped Entrez IDs:\n")
print(length(entrez_ids))

############################################################
# Convert Entrez IDs to gene symbols / names
############################################################

gene_table <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = as.character(entrez_ids),
  keytype = "ENTREZID",
  columns = c("SYMBOL", "GENENAME")
)

gene_table <- unique(gene_table)
gene_table <- gene_table[!is.na(gene_table$SYMBOL), ]

cat("\nMapped gene table preview:\n")
print(head(gene_table))

write.csv(
  gene_table,
  file.path(OUT_DIR, paste0(OUT_PREFIX, "_mapped_genes.csv")),
  row.names = FALSE
)

write.csv(
  unique(gene_table["SYMBOL"]),
  file.path(OUT_DIR, paste0(OUT_PREFIX, "_unique_gene_symbols.csv")),
  row.names = FALSE
)

############################################################
# KEGG enrichment analysis
############################################################

kegg_results <- gometh(
  sig.cpg = top_cpgs,
  all.cpg = all_cpgs,
  collection = "KEGG",
  prior.prob = TRUE,
  array.type = "EPIC"
)

kegg_results <- kegg_results[!is.na(kegg_results$FDR), ]
kegg_results <- kegg_results[is.finite(kegg_results$FDR), ]
kegg_results <- kegg_results[kegg_results$FDR > 0, ]
kegg_results <- kegg_results[order(kegg_results$FDR, -kegg_results$DE), ]

cat("\nTop KEGG pathway results:\n")
print(head(kegg_results, 20))

cat("\nSummary of KEGG FDR values:\n")
print(summary(kegg_results$FDR))

cat("\nSummary of KEGG raw p-values (P.DE):\n")
print(summary(kegg_results$P.DE))

write.csv(
  kegg_results,
  file.path(OUT_DIR, paste0(OUT_PREFIX, "_KEGG_enrichment_all.csv")),
  row.names = FALSE
)

############################################################
# Significant KEGG terms only
############################################################

sig_kegg <- kegg_results[kegg_results$FDR < 0.05, ]

write.csv(
  sig_kegg,
  file.path(OUT_DIR, paste0(OUT_PREFIX, "_KEGG_enrichment_FDRlt0.05.csv")),
  row.names = FALSE
)

cat("\nNumber of significant KEGG pathways (FDR < 0.05):\n")
print(nrow(sig_kegg))

############################################################
# Plot top KEGG pathways
############################################################
kegg_plot_df <- as.data.frame(kegg_results, stringsAsFactors = FALSE)

kegg_plot_df$TERM <- as.character(kegg_plot_df$Description)
kegg_plot_df$N <- suppressWarnings(as.numeric(kegg_plot_df$N))
kegg_plot_df$DE <- suppressWarnings(as.numeric(kegg_plot_df$DE))
kegg_plot_df$P.DE <- suppressWarnings(as.numeric(kegg_plot_df$P.DE))
kegg_plot_df$FDR <- suppressWarnings(as.numeric(kegg_plot_df$FDR))

kegg_plot_df <- kegg_plot_df[!is.na(kegg_plot_df$TERM) & kegg_plot_df$TERM != "", ]
kegg_plot_df <- kegg_plot_df[!is.na(kegg_plot_df$DE) & is.finite(kegg_plot_df$DE), ]
kegg_plot_df <- kegg_plot_df[!is.na(kegg_plot_df$P.DE) & is.finite(kegg_plot_df$P.DE), ]
kegg_plot_df <- kegg_plot_df[!is.na(kegg_plot_df$FDR) & is.finite(kegg_plot_df$FDR), ]
kegg_plot_df <- kegg_plot_df[kegg_plot_df$FDR > 0, ]
kegg_plot_df <- kegg_plot_df[kegg_plot_df$P.DE > 0, ]

kegg_plot_df$qscore <- -log10(kegg_plot_df$FDR)
kegg_plot_df$pscore <- -log10(kegg_plot_df$P.DE)

kegg_plot_df_filtered <- kegg_plot_df[kegg_plot_df$N < 1000, ]

if (nrow(kegg_plot_df_filtered) == 0) {
  warning("No KEGG pathways left after filtering on N < 1000, using unfiltered kegg_plot_df instead.")
  kegg_plot_df_filtered <- kegg_plot_df
}

kegg_plot_df_filtered <- kegg_plot_df_filtered[order(kegg_plot_df_filtered$FDR, -kegg_plot_df_filtered$DE), ]
kegg_plot_df_filtered <- head(kegg_plot_df_filtered, 15)
kegg_plot_df_filtered$TERM_WRAPPED <- str_wrap(kegg_plot_df_filtered$TERM, width = 45)

kegg_plot_df_count <- kegg_plot_df_filtered[order(kegg_plot_df_filtered$DE), ]
kegg_plot_df_count$TERM_WRAPPED <- factor(
  kegg_plot_df_count$TERM_WRAPPED,
  levels = kegg_plot_df_count$TERM_WRAPPED
)

kegg_plot_df_qscore <- kegg_plot_df_filtered[order(kegg_plot_df_filtered$qscore), ]
kegg_plot_df_qscore$TERM_WRAPPED <- factor(
  kegg_plot_df_qscore$TERM_WRAPPED,
  levels = kegg_plot_df_qscore$TERM_WRAPPED
)

kegg_plot_df_pscore <- kegg_plot_df_filtered[order(kegg_plot_df_filtered$pscore), ]
kegg_plot_df_pscore$TERM_WRAPPED <- factor(
  kegg_plot_df_pscore$TERM_WRAPPED,
  levels = kegg_plot_df_pscore$TERM_WRAPPED
)

cat("\nKEGG plot dataframe preview:\n")
print(kegg_plot_df_filtered[, c("TERM", "N", "DE", "P.DE", "FDR", "pscore", "qscore")])
############################################################
# Plot 1: by Count
############################################################

p_kegg_count <- ggplot(
  kegg_plot_df_count,
  aes(x = DE, y = TERM_WRAPPED, fill = qscore)
) +
  geom_col(width = 0.8) +
  labs(
    x = "Number of significant CpGs in pathway",
    y = NULL,
    fill = "-log10(FDR)",
    title = "Top KEGG pathways by count"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 10)
  )

print(p_kegg_count)

ggsave(
  filename = file.path(OUT_DIR, paste0(OUT_PREFIX, "_KEGG_barplot_count.png")),
  plot = p_kegg_count,
  width = 11,
  height = 7,
  dpi = 300
)

############################################################
# Plot 2: by FDR qscore
############################################################

p_kegg_qscore <- ggplot(
  kegg_plot_df_qscore,
  aes(x = qscore, y = TERM_WRAPPED, fill = qscore)
) +
  geom_col(width = 0.8) +
  labs(
    x = expression(-log[10](FDR)),
    y = NULL,
    fill = "-log10(FDR)",
    title = "Top KEGG pathways by FDR significance"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 10)
  )

print(p_kegg_qscore)

ggsave(
  filename = file.path(OUT_DIR, paste0(OUT_PREFIX, "_KEGG_barplot_qscore.png")),
  plot = p_kegg_qscore,
  width = 11,
  height = 7,
  dpi = 300
)

############################################################
# Plot 3: by raw p-value
############################################################

p_kegg_pscore <- ggplot(
  kegg_plot_df_pscore,
  aes(x = pscore, y = TERM_WRAPPED, fill = pscore)
) +
  geom_col(width = 0.8) +
  labs(
    x = expression(-log[10](P.DE)),
    y = NULL,
    fill = "-log10(P.DE)",
    title = "Top KEGG pathways by raw p-value"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 10)
  )

print(p_kegg_pscore)

ggsave(
  filename = file.path(OUT_DIR, paste0(OUT_PREFIX, "_KEGG_barplot_pscore.png")),
  plot = p_kegg_pscore,
  width = 11,
  height = 7,
  dpi = 300
)

############################################################
# Extra checks
############################################################

length(top_cpgs)
length(all_cpgs)
setequal(top_cpgs, all_cpgs)
summary(kegg_results$FDR)
summary(kegg_results$P.DE)

length(unique(gene_table$SYMBOL))
length(entrez_ids)
length(unique(top_cpgs)) / length(unique(gene_table$SYMBOL))

cpg2gene <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = as.character(entrez_ids),
  keytype = "ENTREZID",
  columns = c("SYMBOL")
)

gene_counts <- table(gene_table$SYMBOL)

summary(gene_counts)




head(kegg_results)
colnames(kegg_results)
