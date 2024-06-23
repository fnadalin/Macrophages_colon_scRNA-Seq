
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
	cat("\nUsage: Rscript update_seurat_object.R <in_object> <out_object>\n")
	q()
}

IN_OBJECT <- args[1]
OUT_OBJECT <- args[2]

library("Seurat")

l <- load(IN_OBJECT)
object <- eval(parse(text=l))
object <- UpdateSeuratObject(object)

save(object, file = OUT_OBJECT)

q()

