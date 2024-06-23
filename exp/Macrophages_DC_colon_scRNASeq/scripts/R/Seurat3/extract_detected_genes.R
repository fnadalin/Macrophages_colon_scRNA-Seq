# Created on 01/27/2020

# For each cluster of cells, extract the genes that are detected in >= NUM_CELLS cells in each group

NUM_CELLS <- 20
PERC_CELLS <- 0.1

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
	cat("\nUsage: <obj.Robj> <out_prefix> <meta.slot> [<num_cells> <perc_cells>]\n")
	cat("\n<obj.Robj>   Seurat object created with Seurat v3\n")
	cat("<out_prefix> output file prefix to store detected genes\n")
	cat("<meta.slot>  slot in object@meta.data from where to extract the cell groups\n")
	cat("<num_cells>  minimum number of cells to call a gene as detected in a group\n")
	cat("<perc_cells> minimum fraction of cells to call a gene as detected in a group (in [0,1])\n\n")
	q()
}

LIB_PATH <- file.path(.libPaths()[1], "Seurat_3.0.0")
library("Seurat", lib.loc = LIB_PATH)
library("Matrix")

OBJECT <- args[1]
OUT_PREFIX <- args[2]
META_SLOT <- args[3]
if (length(args) > 3) {
	NUM_CELLS <- as.numeric(args[4])
	if (length(args) > 4) {
		PERC_CELLS <- as.numeric(args[5])
	}
}

object <- eval(parse(text = load(OBJECT)))
M <- GetAssayData(object, slot = "counts")

meta <- object@meta.data[[META_SLOT]]
groups <- unique(meta)

for (g in groups) {
	g_size <- sum(meta == g)
	v <- rowSums(M[,meta == g] > 0)
	names <- rownames(object)[v >= NUM_CELLS & v >= PERC_CELLS*g_size]
	file_name <- paste(OUT_PREFIX, g, ".txt", sep = "")
	write.table(names, file = file_name, quote = FALSE, col.names = FALSE, row.names = FALSE)
}

q()

