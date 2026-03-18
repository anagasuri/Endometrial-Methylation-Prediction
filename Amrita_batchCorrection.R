# subset SH-Data Annotation using R01_Study_IDs_to_include_220622.txt 
library(readxl)

M <- readRDS("M_batch1_batch2.rds")          # CpGs in col1, sample barcodes in other cols

keep_ids <- readLines("R01_Study_IDs_to_include_220622.txt")

# Remove header if present
keep_ids <- keep_ids[keep_ids != "Study.ID"]

annot <- read_excel("SH-Data Annotation_07.07.20_1.xlsx", sheet = "Batch 1and2", col_types = "text")
sum(is.na(annot$`Study ID`))
head(annot$`Study ID`, 20)
tail(annot$`Study ID`, 20)


keep_ids <- readLines("R01_Study_IDs_to_include_220622.txt")
keep_ids <- trimws(keep_ids)
keep_ids <- keep_ids[keep_ids != "Study.ID"]

annot$`Study ID` <- trimws(annot$`Study ID`)

annot_sub <- annot[annot$`Study ID` %in% keep_ids, ]
nrow(annot_sub)


# find epic barcodes from annot_sub
barcodes <- trimws(annot_sub$`Epic_Complete Bar code`)
length(barcodes)          # should be 984
head(barcodes)

# subset methylation matrix M by those barcodes
M_sub <- M[, c(1, match(barcodes, colnames(M))), drop = FALSE]
dim(M_sub)   # should be ~760k x (984 + 1)

# convert to numeric matrix for PCA and SVA
cpg_ids <- M_sub[[1]]
M_mat <- as.matrix(M_sub[, -1, drop = FALSE])
rownames(M_mat) <- cpg_ids

dim(M_mat)   # ~760k x 984

# free memory 
rm(M, M_sub)
gc()

# make sure annotations file aligns with matrix cols
covariates <- as.data.frame(annot_sub)
rownames(covariates) <- barcodes
covariates <- covariates[colnames(M_mat), , drop = FALSE]


stopifnot(identical(rownames(covariates), colnames(M_mat)))

# PCA pre-SVA
pca_pre <- prcomp(t(M_mat), center = FALSE, scale. = FALSE)
dim(pca_pre$x)
summary(pca_pre)
object.size(pca_pre)

# write the PCs to a file 
write.csv(pca_pre$x, "pca_scores.csv")

scores <- read.csv("pca_scores.csv")
head(scores) 

# PCA Plots
inst <- as.factor(covariates[["Institute for Analysis"]])
cols <- rainbow(length(levels(inst)))[inst]

plot(pca_pre$x[,1], pca_pre$x[,2],
     pch=16, col=cols,
     xlab="PC1", ylab="PC2",
     main="Pre-SVA PCA: Institute for Analysis")

legend("topright", legend=levels(inst),
       col=rainbow(length(levels(inst))), pch=16, cex=0.6)

orig <- as.factor(covariates[["Sample Originating Institute"]])
cols <- rainbow(length(levels(orig)))[orig]

plot(pca_pre$x[,1], pca_pre$x[,2],
     pch=16, col=cols,
     xlab="PC1", ylab="PC2",
     main="Pre-SVA PCA: Sample Originating Institute")

legend("topright", legend=levels(orig),
       col=rainbow(length(levels(orig))), pch=16, cex=0.6)

# save plots
dir.create("/home/ubuntu/endo-me_data/pre-SVA_PCA_plots", showWarnings = FALSE)

png("/home/ubuntu/endo-me_data/pre-SVA_PCA_plots/PCA_pre_SVA_Institute_for_Analysis.png",
    width = 2000, height = 1600, res = 300)

inst <- as.factor(covariates[["Institute for Analysis"]])
cols <- rainbow(length(levels(inst)))[inst]

plot(pca_pre$x[,1], pca_pre$x[,2],
     pch=16, col=cols,
     xlab="PC1", ylab="PC2",
     main="Pre-SVA PCA: Institute for Analysis")

legend("topright",
       legend=levels(inst),
       col=rainbow(length(levels(inst))),
       pch=16,
       cex=0.6)

dev.off()

png("/home/ubuntu/endo-me_data/pre-SVA_PCA_plots/PCA_pre_SVA_Batch.png",
    width = 2000, height = 1600, res = 300)

batch <- as.factor(covariates[["Batch"]])
cols <- rainbow(length(levels(batch)))[batch]

plot(pca_pre$x[,1], pca_pre$x[,2],
     pch=16, col=cols,
     xlab="PC1", ylab="PC2",
     main="Pre-SVA PCA: Batch")

legend("topright",
       legend=levels(batch),
       col=rainbow(length(levels(batch))),
       pch=16,
       cex=0.6)

dev.off()

list.files("/home/ubuntu/endo-me_data/pre-SVA_PCA_plots")


# SVA 
stopifnot(identical(rownames(covariates), colnames(M_mat)))

rownames(covariates) <- covariates[["Epic_Complete Bar code"]]
covariates <- covariates[colnames(M_mat), , drop = FALSE]
stopifnot(identical(rownames(covariates), colnames(M_mat)))

covariates$`Endometriosis (Yes/No)` <- factor(covariates$`Endometriosis (Yes/No)`)
covariates$`Cycle phase for Analysis` <- factor(covariates$`Cycle phase for Analysis`)

# install.packages("BiocManager")
# BiocManager::install("sva")        # for EstDimRMT (depending on your setup)
# install.packages("SmartSVA")  # if available in your environment / or source it as you did before

library(sva)

mod   <- model.matrix(~ `Endometriosis (Yes/No)` + `Cycle phase for Analysis`, data = covariates)
mod0  <- model.matrix(~ `Cycle phase for Analysis`, data = covariates)

n.sv <- num.sv(M_mat, mod, method = "leek")
n.sv
# 2 surrogate variables found 

# run SVA and save to rds file 
svobj <- sva(M_mat, mod, mod0, n.sv = n.sv)
saveRDS(svobj, file = "SVA_svobj.rds")

head(svobj)
getwd()
dim(svobj$sv)
# 984 x 2 

# 

#########################################
#########################################
#########################################

############################################################
# Endo methylation: subset -> PCA pre -> SmartSVA -> correct -> PCA post
# Saves: same pre-SVA plots + PCA score CSVs + SmartSVA objects + corrected matrix + post-SVA plots
############################################################

# subset SH-Data Annotation using R01_Study_IDs_to_include_220622.txt 
suppressPackageStartupMessages({
  library(readxl)
  library(sva)
  library(SmartSVA)
})

# -----------------------------
# 0) Load methylation matrix
# -----------------------------
M <- readRDS("M_batch1_batch2.rds")  # CpGs in col1, sample barcodes in other cols

# -----------------------------
# 1) Read keep IDs + annotation, subset to study IDs
# -----------------------------
keep_ids <- readLines("R01_Study_IDs_to_include_220622.txt")
keep_ids <- trimws(keep_ids)
keep_ids <- keep_ids[keep_ids != "Study.ID"]  # remove header if present

annot <- read_excel("SH-Data Annotation_07.07.20_1.xlsx",
                    sheet = "Batch 1and2",
                    col_types = "text")

annot$`Study ID` <- trimws(annot$`Study ID`)
annot_sub <- annot[annot$`Study ID` %in% keep_ids, ]
nrow(annot_sub)

# -----------------------------
# 2) Find EPIC barcodes and subset methylation matrix
# -----------------------------
barcodes <- trimws(annot_sub$`Epic_Complete Bar code`)
length(barcodes)           # should be 984
head(barcodes)

# Keep CpG ID column + barcode columns (in the barcode order)
idx <- match(barcodes, colnames(M))

# Fail early if any barcodes are missing in M
if (any(is.na(idx))) {
  missing_barcodes <- barcodes[is.na(idx)]
  stop("These barcodes are in annotation but missing from M colnames:\n",
       paste(missing_barcodes, collapse = "\n"))
}

M_sub <- M[, c(1, idx), drop = FALSE]
dim(M_sub)                 # ~760k x (984 + 1)

# -----------------------------
# 3) Convert to numeric matrix (CpGs x samples)
# -----------------------------
cpg_ids <- M_sub[[1]]
M_mat <- as.matrix(M_sub[, -1, drop = FALSE])
rownames(M_mat) <- cpg_ids
dim(M_mat)                 # ~760k x 984

# free memory
rm(M, M_sub)
gc() 

# -----------------------------
# 4) Align covariates rows to M_mat columns
# -----------------------------
covariates <- as.data.frame(annot_sub)
rownames(covariates) <- barcodes
covariates <- covariates[colnames(M_mat), , drop = FALSE]
stopifnot(identical(rownames(covariates), colnames(M_mat)))

# -----------------------------
# 5) PCA pre-SVA (same as before) + save outputs
# -----------------------------
pca_pre <- prcomp(t(M_mat), center = FALSE, scale. = FALSE)

# write PCs to the same file name as before
write.csv(pca_pre$x, "pca_scores.csv")

scores <- read.csv("pca_scores.csv")
head(scores)

# Plot to screen (same as before)
inst <- as.factor(covariates[["Institute for Analysis"]])
cols <- rainbow(length(levels(inst)))[inst]
plot(pca_pre$x[,1], pca_pre$x[,2],
     pch=16, col=cols,
     xlab="PC1", ylab="PC2",
     main="Pre-SVA PCA: Institute for Analysis")
legend("topright", legend=levels(inst),
       col=rainbow(length(levels(inst))), pch=16, cex=0.6)

orig <- as.factor(covariates[["Sample Originating Institute"]])
cols <- rainbow(length(levels(orig)))[orig]
plot(pca_pre$x[,1], pca_pre$x[,2],
     pch=16, col=cols,
     xlab="PC1", ylab="PC2",
     main="Pre-SVA PCA: Sample Originating Institute")
legend("topright", legend=levels(orig),
       col=rainbow(length(levels(orig))), pch=16, cex=0.6)

# Save plots (same directory + filenames as before)
dir.create("/home/ubuntu/endo-me_data/pre-SVA_PCA_plots", showWarnings = FALSE, recursive = TRUE)

png("/home/ubuntu/endo-me_data/pre-SVA_PCA_plots/PCA_pre_SVA_Institute_for_Analysis.png",
    width = 2000, height = 1600, res = 300)
inst <- as.factor(covariates[["Institute for Analysis"]])
cols <- rainbow(length(levels(inst)))[inst]
plot(pca_pre$x[,1], pca_pre$x[,2],
     pch=16, col=cols,
     xlab="PC1", ylab="PC2",
     main="Pre-SVA PCA: Institute for Analysis")
legend("topright",
       legend=levels(inst),
       col=rainbow(length(levels(inst))),
       pch=16,
       cex=0.6)
dev.off()

png("/home/ubuntu/endo-me_data/pre-SVA_PCA_plots/PCA_pre_SVA_Batch.png",
    width = 2000, height = 1600, res = 300)
batch <- as.factor(covariates[["Batch"]])
cols <- rainbow(length(levels(batch)))[batch]
plot(pca_pre$x[,1], pca_pre$x[,2],
     pch=16, col=cols,
     xlab="PC1", ylab="PC2",
     main="Pre-SVA PCA: Batch")
legend("topright",
       legend=levels(batch),
       col=rainbow(length(levels(batch))),
       pch=16,
       cex=0.6)
dev.off()

list.files("/home/ubuntu/endo-me_data/pre-SVA_PCA_plots")

# -----------------------------
# 6) SmartSVA (Idit-style) + save objects
# -----------------------------
# Make sure rownames match sample IDs (same checks as before)
stopifnot(identical(rownames(covariates), colnames(M_mat)))

# Ensure factors
covariates$`Endometriosis (Yes/No)` <- factor(covariates$`Endometriosis (Yes/No)`)
covariates$`Cycle phase for Analysis` <- factor(covariates$`Cycle phase for Analysis`)

# Models: match analysis_model_SVA_final.R (mod keeps biology; mod0 is intercept-only)
mod  <- model.matrix(~ `Endometriosis (Yes/No)` + `Cycle phase for Analysis`, data = covariates)
mod0 <- model.matrix(~ 1, data = covariates)

# Estimate number of SVs like the shared script (EstDimRMT on residualized data)
Y.r <- t(resid(lm(t(M_mat) ~ `Endometriosis (Yes/No)` + `Cycle phase for Analysis`, data = covariates)))
n.sv <- EstDimRMT(Y.r, FALSE)$dim + 1
n.sv

# Run SmartSVA
svaout <- smartsva.cpp(as.matrix(M_mat), n.sv = n.sv, mod = mod, mod0 = mod0)
svaout_sv <- svaout$sv
rownames(svaout_sv) <- colnames(M_mat)
dim(svaout_sv)

# Save SmartSVA outputs (RDS in current working directory, like you did for SVA_svobj.rds)
saveRDS(svaout, file = "SmartSVA_svaout.rds")
saveRDS(svaout_sv, file = "SmartSVA_sv.rds")

# -----------------------------
# 7) Correct matrix (exactly like the shared script) + save
# -----------------------------
get_corrected_data <- function(y, mod, svaobj) {
  X <- cbind(mod, svaobj$sv)
  Hat <- solve(t(X) %*% X) %*% t(X)
  beta <- (Hat %*% t(y))
  P <- ncol(mod)
  cleany <- y - t(as.matrix(X[, -(1:P), drop = FALSE]) %*% beta[-(1:P), , drop = FALSE])
  return(cleany)
}

M_corrected <- get_corrected_data(as.matrix(M_mat), mod, svaout)

# Save corrected matrix (RDS in current working directory)
saveRDS(M_corrected, file = "M_mat_SmartSVA_corrected.rds")

# -----------------------------
# 8) PCA post-correction + save to CSV + plots
# -----------------------------
pca_post <- prcomp(t(M_corrected), center = FALSE, scale. = FALSE)

# Save post PCs (new file; pre file name preserved as pca_scores.csv)
write.csv(pca_post$x, "pca_scores_post_SmartSVA.csv")

# Save post plots (new directory; pre directory preserved)
dir.create("/home/ubuntu/endo-me_data/post-SVA_PCA_plots", showWarnings = FALSE, recursive = TRUE)

png("/home/ubuntu/endo-me_data/post-SVA_PCA_plots/PCA_post_SVA_Institute_for_Analysis.png",
    width = 2000, height = 1600, res = 300)
inst <- as.factor(covariates[["Institute for Analysis"]])
cols <- rainbow(length(levels(inst)))[inst]
plot(pca_post$x[,1], pca_post$x[,2],
     pch=16, col=cols,
     xlab="PC1", ylab="PC2",
     main="Post-SmartSVA PCA: Institute for Analysis")
legend("topright",
       legend=levels(inst),
       col=rainbow(length(levels(inst))),
       pch=16,
       cex=0.6)
dev.off()

png("/home/ubuntu/endo-me_data/post-SVA_PCA_plots/PCA_post_SVA_Batch.png",
    width = 2000, height = 1600, res = 300)
batch <- as.factor(covariates[["Batch"]])
cols <- rainbow(length(levels(batch)))[batch]
plot(pca_post$x[,1], pca_post$x[,2],
     pch=16, col=cols,
     xlab="PC1", ylab="PC2",
     main="Post-SmartSVA PCA: Batch")
legend("topright",
       legend=levels(batch),
       col=rainbow(length(levels(batch))),
       pch=16,
       cex=0.6)
dev.off()

list.files("/home/ubuntu/endo-me_data/post-SVA_PCA_plots")

# -----------------------------
# 9) Quick checks
# -----------------------------
print(getwd())
print(dim(M_mat))
print(dim(M_corrected))
print(dim(svaout_sv))



head(SmartSVA_svaout)
head(SmartSVA_sv)
