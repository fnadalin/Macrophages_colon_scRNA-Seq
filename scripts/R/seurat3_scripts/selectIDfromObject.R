# Created on 2019/06/17

# Select the raw and normalized data from a list of IDs and create a new object

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
	cat("\nUsage: Rscript selectIDfromObject.R <in_object> <identity_slot> <IDs> <out_object>\n")
	cat("<out_object>         input Seurat object file name\n")
	cat("<identity_slot>      name of the identity slot in the object\n")
	cat("<IDs>                comma-separated list of the IDs to select\n")
	cat("<out_object>         output Seurat object file name\n\n")
	q()
}

IN_OBJECT <- args[1]
IDENTITY_SLOT <- args[2]
IDs <- args[3]
OUT_OBJECT <- args[4]

library("Seurat")

object <- eval(parse(text=load(IN_OBJECT)))
ids <- unlist(strsplit(IDs, split=","))
object <- SetIdent(object, value = IDENTITY_SLOT)
object <- SubsetData(object, ident.use = ids, do.clean = FALSE, subset.raw = TRUE, do.center = TRUE, do.scale = TRUE)
object <- SetIdent(object, value = "orig.ident")
save(object, file = OUT_OBJECT)

q()

