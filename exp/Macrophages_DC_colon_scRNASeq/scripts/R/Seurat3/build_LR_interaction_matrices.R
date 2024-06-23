# Created on 01/27/2020

# Collect the results of count_LR_interactions.R in a single matrix of absolute counts (LR and RL)
# Try also some normalization formulas...

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
	cat("\nUsage: <dir> <num_cl_1> <num_cl_2> <out_prefix>\n")
	cat("\n<dir>        directory where the results of count_LR_interactions.R are stored\n")
	cat("<num_cl_1>     number of clusters for set1 (MP)\n")
	cat("<num_cl_2>     number of clusters for set2 (DC)\n")
	cat("<out_prefix>   for the matrices (abs and norm)\n\n")
	q()
}

DIR <- args[1]
NUM_CL_1 <- as.numeric(args[2])
NUM_CL_2 <- as.numeric(args[3])
OUT_PREFIX <- args[4]

M_LR <- matrix(0, nrow = NUM_CL_1, ncol = NUM_CL_2)
M_RL <- matrix(0, nrow = NUM_CL_1, ncol = NUM_CL_2)
M_LR_norm1 <- matrix(0, nrow = NUM_CL_1, ncol = NUM_CL_2)
M_RL_norm1 <- matrix(0, nrow = NUM_CL_1, ncol = NUM_CL_2)

M_LR <- apply(M_LR, 2, as.numeric)
M_RL <- apply(M_RL, 2, as.numeric)

for (i in 0:(NUM_CL_1-1)) {
	for (j in 0:(NUM_CL_2-1)) {
		file_name <- file.path(DIR, paste("interactions_cl", i, "-cl", j, ".txt", sep=""))
		df <- as.matrix(read.table(file_name, sep="="))
		num_genes1_L <- as.numeric(df[df[,1] == "num_genes1_L",2])
		num_genes1_R <- as.numeric(df[df[,1] == "num_genes1_R",2])
		num_genes2_L <- as.numeric(df[df[,1] == "num_genes2_L",2])
		num_genes2_R <- as.numeric(df[df[,1] == "num_genes2_R",2])
		M_LR[i+1,j+1] <- M_LR_norm1[i+1,j+1] <- as.numeric(df[df[,1] == "LR_pairs",2])
		M_RL[i+1,j+1] <- M_RL_norm1[i+1,j+1] <- as.numeric(df[df[,1] == "RL_pairs",2])
		if (M_LR[i+1,j+1] > 0) {
			M_LR_norm1[i+1,j+1] <- M_LR[i+1,j+1]/(num_genes1_L*num_genes2_R)
		} 
		if (M_RL[i+1,j+1] > 0) {
			M_RL_norm1[i+1,j+1] <- M_RL[i+1,j+1]/(num_genes1_R*num_genes2_L)
		}
	}
}

rownames(M_LR) <- rownames(M_LR_norm1) <- as.character(0:(NUM_CL_1-1))
rownames(M_RL) <- rownames(M_RL_norm1) <- as.character(0:(NUM_CL_1-1))
colnames(M_LR) <- colnames(M_LR_norm1) <- as.character(0:(NUM_CL_2-1))
colnames(M_RL) <- colnames(M_RL_norm1) <- as.character(0:(NUM_CL_2-1))

write.table(M_LR, file = paste(OUT_PREFIX, "_LR.tsv", sep=""), quote = FALSE, sep = "\t")
write.table(M_RL, file = paste(OUT_PREFIX, "_RL.tsv", sep=""), quote = FALSE, sep = "\t")
write.table(M_LR_norm1, file = paste(OUT_PREFIX, "_LR_norm1.tsv", sep=""), quote = FALSE, sep = "\t")
write.table(M_RL_norm1, file = paste(OUT_PREFIX, "_RL_norm1.tsv", sep=""), quote = FALSE, sep = "\t")


q()



