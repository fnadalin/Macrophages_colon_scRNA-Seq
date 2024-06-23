# Created on 01/08/2020

# Convert an RDS object into a saved object (Seurat 3)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
	cat("\nUsage: <obj.rds> <obj.Robj>\n")
	q()
}

OBJECT_RDS <- args[1]
OBJECT_ROBJ <- args[2]

library("Seurat")

object <- readRDS(OBJECT_RDS)
save(object, file = OBJECT_ROBJ)

q()

