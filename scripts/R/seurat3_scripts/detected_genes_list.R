
# Created on 01/20/2020

# count the number of genes detected in a threshold number of cells

ORG <- "Hs"

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
	cat("\nUsage: Rscript detected_genes.R <object.Robj> <perc_cells> <max_cells> <out_file> [<org_id>]\n\n")
	cat("<org_id> = \"Hs\" (default) or \"Mm\"\n\n")
	q()
}

OBJECT <- args[1]
PERC_CELLS <- as.numeric(args[2])
MAX_CELLS <- as.numeric(args[3])
OUT_FILE <- args[4]
if (length(args) > 4) {
	ORG <- args[5]
}

library("Seurat")
library("Matrix")
library("limma")

object <- eval(parse(text=load(OBJECT)))
M <- GetAssayData(object, slot = "counts")
NUM_CELLS <- floor(min(MAX_CELLS,ncol(M)*PERC_CELLS/100))
cat(paste("Select the genes detected in >=", NUM_CELLS, "cells\n"))

genes <- rownames(M)[rowSums(M != 0) >= NUM_CELLS]
genes_official <- alias2SymbolTable(genes, species=ORG)
df <- data.frame(genes, genes_official)
write.table(df, file=OUT_FILE, row.names=FALSE, sep="\t", quote=FALSE)

q()

