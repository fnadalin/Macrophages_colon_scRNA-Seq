
# Created on 02/07/2020

# count the number of genes detected in a threshold number of cells

IDENT_SLOT <- "orig.ident"
IDENTS <- NULL

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
	cat("\nUsage: Rscript detected_genes.R <object.Robj> <num_cells> [<ident_slot> <idents>]\n\n")
	cat("<ident_slot>     identity slot to extract groups of cells from\n")
	cat("<idents>         comma-separated identities identifying the groups of cells to extract\n\n")
	q()
}

OBJECT <- args[1]
NUM_CELLS <- as.numeric(args[2])
if (length(args) > 2) {
	IDENT_SLOT <- args[3]
	if (length(args) > 3) {
		IDENTS <- unlist(strsplit(args[4], split=","))
	}
}

library("Seurat")
library("Matrix")

object <- eval(parse(text=load(OBJECT)))
M <- GetAssayData(object, slot = "counts")
if (!is.null(IDENTS)) {
	idx <- object@meta.data[[IDENT_SLOT]] %in% IDENTS
	M <- M[,idx]
}
genes <- rownames(M)[rowSums(M != 0) >= NUM_CELLS]
cat(paste(length(genes), "\n", sep=""))

q()

