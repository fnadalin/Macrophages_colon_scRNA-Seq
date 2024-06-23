# Created on 01/28/2020

# Plot the interaction matrices as annotated heatmaps

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
	cat("\nUsage: <out_prefix>\n")
	cat("\n<dir>        directory where the results of count_LR_interactions.R are stored\n")
	cat("<out_prefix>   for the matrices (abs and norm)\n\n")
	q()
}

OUT_PREFIX <- args[1]

library("ComplexHeatmap")
library("scales")

M_LR <- as.matrix(read.table(paste(OUT_PREFIX, "_LR.tsv", sep=""), sep = "\t", header = TRUE))
M_RL <- as.matrix(read.table(paste(OUT_PREFIX, "_RL.tsv", sep=""), sep = "\t", header = TRUE))
M_LR_norm1 <- as.matrix(read.table(paste(OUT_PREFIX, "_LR_norm1.tsv", sep=""), sep = "\t", header = TRUE))
M_RL_norm1 <- as.matrix(read.table(paste(OUT_PREFIX, "_RL_norm1.tsv", sep=""), sep = "\t", header = TRUE))

colnames(M_LR) <- colnames(M_RL) <- colnames(M_LR_norm1) <- colnames(M_RL_norm1) <- gsub("^X","",colnames(M_LR))

ht_LR <- Heatmap(M_LR, name = "LR count", show_row_names = TRUE, show_column_names = TRUE, row_title = " ", column_title = "ligand-receptor", 
		cluster_rows = FALSE, cluster_columns = FALSE, column_title_gp = gpar(fontsize = 10), row_names_side = "left", column_names_side = "bottom")
ht_RL <- Heatmap(M_RL, name = "RL count", show_row_names = FALSE, show_column_names = TRUE, row_title = " ", column_title = "receptor-ligand", 
		cluster_rows = FALSE, cluster_columns = FALSE, column_title_gp = gpar(fontsize = 10), row_names_side = "left", column_names_side = "bottom")

ht_LR_norm <- Heatmap(M_LR_norm1, name = "LR norm. count", show_row_names = TRUE, show_column_names = TRUE, row_title = " ", column_title = "ligand-receptor", 
		cluster_rows = FALSE, cluster_columns = FALSE, column_title_gp = gpar(fontsize = 10), row_names_side = "left", column_names_side = "bottom")
ht_RL_norm <- Heatmap(M_RL_norm1, name = "RL norm. count", show_row_names = FALSE, show_column_names = TRUE, row_title = " ", column_title = "receptor-ligand", 
		cluster_rows = FALSE, cluster_columns = FALSE, column_title_gp = gpar(fontsize = 10), row_names_side = "left", column_names_side = "bottom")

pdf(paste(OUT_PREFIX, "_count_heatmap.pdf", sep=""), width=2+0.6*ncol(M_LR), height=0.7*nrow(M_LR))
draw(ht_LR + ht_RL, row_title = "MP", column_title = "DC", row_title_gp = gpar(fontface = "bold", fontsize = 16), column_title_gp = gpar(fontface = "bold", fontsize = 16))
dev.off()

pdf(paste(OUT_PREFIX, "_norm_heatmap.pdf", sep=""), width=2+0.6*ncol(M_LR), height=0.7*nrow(M_LR))
draw(ht_LR_norm + ht_RL_norm, row_title = "MP", column_title = "DC", row_title_gp = gpar(fontface = "bold", fontsize = 16), column_title_gp = gpar(fontface = "bold", fontsize = 16))
dev.off()


q()
