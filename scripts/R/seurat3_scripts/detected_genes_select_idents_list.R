
# Created on 02/10/2020

# count the number of genes detected in a threshold number of cells

IDENT_SLOT <- "orig.ident"
IDENTS <- NULL

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
	cat("\nUsage: Rscript detected_genes.R <object.Robj> <num_cells> <out_file> [<ident_slot> <idents>]\n\n")
	cat("<ident_slot>     identity slot to extract groups of cells from\n")
	cat("<idents>         comma-separated identities identifying the groups of cells to extract\n\n")
	q()
}

OBJECT <- args[1]
NUM_CELLS <- as.numeric(args[2])
OUT_FILE <- args[3]
if (length(args) > 3) {
	IDENT_SLOT <- args[4]
	if (length(args) > 4) {
		IDENTS <- unlist(strsplit(args[5], split=","))
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
write.table(genes, file=OUT_FILE, row.names=FALSE, col.names=FALSE, quote=FALSE)


q()

