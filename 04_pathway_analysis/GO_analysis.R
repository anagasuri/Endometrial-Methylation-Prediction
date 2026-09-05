############################################################
# CpG -> Entrez gene mapping + GO pathway analysis
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
OUT_DIR <- "GO_analysis"

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
# GO enrichment analysis
############################################################

go_results <- gometh(
  sig.cpg = top_cpgs,
  all.cpg = all_cpgs,
  collection = "GO",
  prior.prob = TRUE,
  array.type = "EPIC"
)

go_results_bp <- go_results[go_results$ONTOLOGY == "BP", ]
go_results_bp <- go_results_bp[!is.na(go_results_bp$FDR), ]
go_results_bp <- go_results_bp[is.finite(go_results_bp$FDR), ]
go_results_bp <- go_results_bp[go_results_bp$FDR > 0, ]
go_results_bp <- go_results_bp[order(go_results_bp$FDR, -go_results_bp$DE), ]

cat("\nTop GO Biological Process results:\n")
print(head(go_results_bp, 20))

cat("\nSummary of GO BP FDR values:\n")
print(summary(go_results_bp$FDR))

cat("\nSummary of GO BP raw p-values (P.DE):\n")
print(summary(go_results_bp$P.DE))

write.csv(
  go_results,
  file.path(OUT_DIR, paste0(OUT_PREFIX, "_GO_enrichment_all.csv")),
  row.names = FALSE
)

write.csv(
  go_results_bp,
  file.path(OUT_DIR, paste0(OUT_PREFIX, "_GO_enrichment_BP.csv")),
  row.names = FALSE
)

############################################################
# Significant GO terms only
############################################################

sig_go_bp <- go_results_bp[go_results_bp$FDR < 0.05, ]

write.csv(
  sig_go_bp,
  file.path(OUT_DIR, paste0(OUT_PREFIX, "_GO_enrichment_BP_FDRlt0.05.csv")),
  row.names = FALSE
)

cat("\nNumber of significant BP GO terms (FDR < 0.05):\n")
print(nrow(sig_go_bp))

############################################################
# Plot top GO Biological Process pathways
############################################################

plot_df <- as.data.frame(go_results_bp, stringsAsFactors = FALSE)

plot_df$TERM <- as.character(plot_df$TERM)
plot_df$ONTOLOGY <- as.character(plot_df$ONTOLOGY)
plot_df$N <- suppressWarnings(as.numeric(plot_df$N))
plot_df$DE <- suppressWarnings(as.numeric(plot_df$DE))
plot_df$P.DE <- suppressWarnings(as.numeric(plot_df$P.DE))
plot_df$FDR <- suppressWarnings(as.numeric(plot_df$FDR))

plot_df <- plot_df[!is.na(plot_df$TERM) & plot_df$TERM != "", ]
plot_df <- plot_df[!is.na(plot_df$DE) & is.finite(plot_df$DE), ]
plot_df <- plot_df[!is.na(plot_df$P.DE) & is.finite(plot_df$P.DE), ]
plot_df <- plot_df[!is.na(plot_df$FDR) & is.finite(plot_df$FDR), ]
plot_df <- plot_df[plot_df$FDR > 0, ]
plot_df <- plot_df[plot_df$P.DE > 0, ]

plot_df$qscore <- -log10(plot_df$FDR)
plot_df$pscore <- -log10(plot_df$P.DE)

plot_df_filtered <- plot_df[plot_df$N < 1000, ]

if (nrow(plot_df_filtered) == 0) {
  warning("No terms left after filtering on N < 1000, using unfiltered plot_df instead.")
  plot_df_filtered <- plot_df
}

plot_df_filtered <- plot_df_filtered[order(plot_df_filtered$FDR, -plot_df_filtered$DE), ]
plot_df_filtered <- head(plot_df_filtered, 15)
plot_df_filtered$TERM_WRAPPED <- str_wrap(plot_df_filtered$TERM, width = 45)

plot_df_count <- plot_df_filtered[order(plot_df_filtered$DE), ]
plot_df_count$TERM_WRAPPED <- factor(
  plot_df_count$TERM_WRAPPED,
  levels = plot_df_count$TERM_WRAPPED
)

plot_df_qscore <- plot_df_filtered[order(plot_df_filtered$qscore), ]
plot_df_qscore$TERM_WRAPPED <- factor(
  plot_df_qscore$TERM_WRAPPED,
  levels = plot_df_qscore$TERM_WRAPPED
)

plot_df_pscore <- plot_df_filtered[order(plot_df_filtered$pscore), ]
plot_df_pscore$TERM_WRAPPED <- factor(
  plot_df_pscore$TERM_WRAPPED,
  levels = plot_df_pscore$TERM_WRAPPED
)

cat("\nGO plot dataframe preview:\n")
print(plot_df_filtered[, c("TERM", "N", "DE", "P.DE", "FDR", "pscore", "qscore")])

############################################################
# Plot 1: by Count
############################################################

p_count <- ggplot(
  plot_df_count,
  aes(x = DE, y = TERM_WRAPPED, fill = qscore)
) +
  geom_col(width = 0.8) +
  labs(
    x = "Number of significant CpGs in term",
    y = NULL,
    fill = "-log10(FDR)",
    title = "Top GO Biological Process pathways by count"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 10)
  )

print(p_count)

ggsave(
  filename = file.path(OUT_DIR, paste0(OUT_PREFIX, "_GO_BP_barplot_count.png")),
  plot = p_count,
  width = 11,
  height = 7,
  dpi = 300
)

############################################################
# Plot 2: by FDR qscore
############################################################

p_qscore <- ggplot(
  plot_df_qscore,
  aes(x = qscore, y = TERM_WRAPPED, fill = qscore)
) +
  geom_col(width = 0.8) +
  labs(
    x = expression(-log[10](FDR)),
    y = NULL,
    fill = "-log10(FDR)",
    title = "Top GO Biological Process pathways by FDR significance"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 10)
  )

print(p_qscore)

ggsave(
  filename = file.path(OUT_DIR, paste0(OUT_PREFIX, "_GO_BP_barplot_qscore.png")),
  plot = p_qscore,
  width = 11,
  height = 7,
  dpi = 300
)

############################################################
# Plot 3: by raw p-value
############################################################

p_pscore <- ggplot(
  plot_df_pscore,
  aes(x = pscore, y = TERM_WRAPPED, fill = pscore)
) +
  geom_col(width = 0.8) +
  labs(
    x = expression(-log[10](P.DE)),
    y = NULL,
    fill = "-log10(P.DE)",
    title = "Top GO Biological Process pathways by raw p-value"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 10)
  )

print(p_pscore)

ggsave(
  filename = file.path(OUT_DIR, paste0(OUT_PREFIX, "_GO_BP_barplot_pscore.png")),
  plot = p_pscore,
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
summary(go_results_bp$FDR)
summary(go_results_bp$P.DE)

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