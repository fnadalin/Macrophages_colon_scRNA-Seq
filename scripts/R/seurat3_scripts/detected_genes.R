
# Created on 01/16/2020

# count the number of genes detected in a threshold number of cells

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
	cat("\nUsage: Rscript detected_genes.R <object.Robj> <num_cells>\n\n")
	q()
}

OBJECT <- args[1]
NUM_CELLS <- as.numeric(args[2])

library("Seurat")
library("Matrix")

object <- eval(parse(text=load(OBJECT)))
M <- GetAssayData(object, slot = "counts")
genes <- rownames(M)[rowSums(M != 0) >= NUM_CELLS]
cat(paste(length(genes), "\n", sep=""))

q()

