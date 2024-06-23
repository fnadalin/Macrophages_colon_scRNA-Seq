# Created on 2019/06/19

# sample the same fraction from each ID from the active.identity slot

FRAC <- 0.5
SEED <- 123

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
	cat("\nUsage: <object.Robj> [<frac_cells> <seed> <meta.slot>]\n")
	cat("<object.Robj> Seurat3 object\n")
	cat("<frac_cells>  fraction of cells to extract from object\n")
	cat("<seed>        seed for the RNG\n")
	cat("<meta.slot>   name of the meta.data slot where the sampling result is stored (cell labelling: sample1 and sample2)\n\n")
	q()
}

OBJECT <- args[1]
if (length(args) > 1) {
	FRAC <- as.numeric(args[2])
	if (length(args) > 2) {
		SEED <- as.numeric(args[3])
	}
}
if (length(args) > 3) {
	META_SLOT <- args[4]
} else {
	META_SLOT <- paste("sampling", SEED, sep=".")
}

library("Seurat")

object <- eval(parse(text=load(OBJECT)))

cell_ids <- Idents(object)
ids <- unique(cell_ids)

set.seed(SEED)
cells <- 1:ncol(GetAssayData(object))

meta <- rep("sample1", length(cells))
for (id in ids) {
	c <- (cell_ids == id)
	s <- sample(cells[c], size=sum(c)*FRAC)
	meta[s] <- "sample2"
}
meta <- paste(Idents(object), meta, sep="-")
object[[META_SLOT]] <- meta # check that this works!!!

save(object, file=OBJECT)

q()

