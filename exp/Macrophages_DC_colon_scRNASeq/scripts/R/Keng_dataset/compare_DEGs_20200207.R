
# Created on 02/07/2020

library("ggplot2")
library("ggrepel")

ANALYSIS_PATH <- "/storage2/Macrophages_DC_colon_scRNASeq/analysis"
BLP_PATH <- "Seurat/optimal_clustering/test7_Seurat3/Macro"
BLP_DEG_PATH <- "Seurat/optimal_clustering/test7_Seurat3/Macro_optimal_clustering/param_set23/Macro/MAST"
KANG_PATH <- "data/Macro_Kang2019"
OUT_DIR <- file.path(ANALYSIS_PATH, "Seurat/Kang_dataset_20200207")

BLP_MP_DEG_FILE <- file.path(ANALYSIS_PATH, BLP_DEG_PATH, "DEG_MAST_cl0-1.tsv")
KANG_MP_DEG_FILE <- file.path(ANALYSIS_PATH, KANG_PATH, "DEG_cl4-cl6.tsv")
DETECTED_GENES <- file.path(OUT_DIR, "detected_genes_20cells.txt")

BLP_MP_DEG <- read.table(BLP_MP_DEG_FILE, header = TRUE)
KANG_MP_DEG <- read.table(KANG_MP_DEG_FILE, header = TRUE)

COMMON_UP <- file.path(OUT_DIR, "Kang_cl4_UP-Sasha_cl1_UP.txt")
COMMON_DOWN <- file.path(OUT_DIR, "Kang_cl6_UP-Sasha_cl0_UP.txt")

# transform the logFC
BLP_MP_DEG$avg_log2FC <- -BLP_MP_DEG$avg_log2FC
KANG_MP_DEG <- cbind(KANG_MP_DEG, avg_log2FC = KANG_MP_DEG$avg_logFC/log(2))
KANG_MP_DEG$avg_log2FC <- -KANG_MP_DEG$avg_log2FC

# extract the UP- and DOWN- regulated genes 
BLP_MP_UP <- BLP_MP_DEG[BLP_MP_DEG$avg_log2FC > 0.5 & BLP_MP_DEG$p_val_adj < 0.05,]
BLP_MP_DOWN <- BLP_MP_DEG[BLP_MP_DEG$avg_log2FC < -0.5 & BLP_MP_DEG$p_val_adj < 0.05,]
KANG_MP_UP <- KANG_MP_DEG[KANG_MP_DEG$avg_log2FC > 0.5 & KANG_MP_DEG$p_val_adj < 0.05,]
KANG_MP_DOWN <- KANG_MP_DEG[KANG_MP_DEG$avg_log2FC < -0.5 & KANG_MP_DEG$p_val_adj < 0.05,]

common_UP <- as.character(BLP_MP_UP$geneID[BLP_MP_UP$geneID %in% KANG_MP_UP$gene])
common_DOWN <- as.character(BLP_MP_DOWN$geneID[BLP_MP_DOWN$geneID %in% KANG_MP_DOWN$gene])

write.table(common_UP, file = COMMON_UP, quote = FALSE, col.names = FALSE, row.names = FALSE)
write.table(common_DOWN, file = COMMON_DOWN, quote = FALSE, col.names = FALSE, row.names = FALSE)

# create the data frame containing the union of the UP-regulated DEGs
gene <- c(as.character(BLP_MP_UP$geneID), as.character(KANG_MP_UP$gene))
avg_log2FC <- c(BLP_MP_UP$avg_log2FC, KANG_MP_UP$avg_log2FC)
experiment <- c(rep("BLP",nrow(BLP_MP_UP)), rep("KANG",nrow(KANG_MP_UP)))
experiment <- factor(experiment, levels = c("KANG", "BLP"))
group <- gene
group[!(gene %in% common_UP)] <- NA
df_up <- data.frame(gene, avg_log2FC, experiment, group)

# run the Fisher tests

universe <- as.numeric(drop(as.matrix(read.table(DETECTED_GENES))))

a_and_b <- length(common_UP)
a_not_b <- nrow(BLP_MP_UP) - length(common_UP)
b_not_a <- nrow(KANG_MP_UP) - length(common_UP)
not_a_not_b <- universe - nrow(BLP_MP_UP) - nrow(KANG_MP_UP) - length(common_UP)
F <- matrix(c(a_and_b, b_not_a, a_not_b, not_a_not_b), nrow=2, ncol=2)
cat("===== UP IN BLP+MP AND CL4 =====\n")
write.table(F)
fisher.test(x = F, alternative = "greater")

a_and_b <- length(common_DOWN)
a_not_b <- nrow(BLP_MP_DOWN) - length(common_DOWN)
b_not_a <- nrow(KANG_MP_DOWN) - length(common_DOWN)
not_a_not_b <- universe - nrow(BLP_MP_UP) - nrow(KANG_MP_DOWN) - length(common_DOWN)
F <- matrix(c(a_and_b, b_not_a, a_not_b, not_a_not_b), nrow=2, ncol=2)
cat("===== UP IN BLP-MP AND CL6 =====\n")
write.table(F)
fisher.test(x = F, alternative = "greater")







q()


