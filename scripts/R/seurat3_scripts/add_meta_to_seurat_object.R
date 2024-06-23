# Created on 2019/10/09

# Add meta.data value field to Seurat object containing 1=doublet, 0=non-doublet

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
	cat("Usage: Rscript add_meta_to_seurat_object.R <values> <cell_IDs> <meta_field> <seurat_object>\n")
	cat("<values>           list of meta-data values, one per line")
	cat("<cell_IDs>         list of cell IDs, one per line (missing cellIDs are assigned meta.data=0)")
	cat("<meta_field>       to be created inside object@meta.data")
	cat("<seurat_object>    where the Seurat object is stored")
	q()
}

VALUES <- args[1]
CELL_ID <- args[2]
META_FIELD <- args[3]
OBJECT <- args[4]

library("Seurat")

obj_name <- load(OBJECT)
object <- eval(parse(text = obj_name))

df1 <- try(read.table(VALUES))
df2 <- try(read.table(CELL_ID))
if (inherits(df1, 'try-error') | inherits(df2, 'try-error')) { 
	values <- cell_id <- NULL
} else {
	values <- c(as.matrix(df1))
	cell_id <- gsub("-1", "", c(as.matrix(df2)))
}

if (length(values) != length(cell_id)) {
	write(paste("length(values)=", length(values), ", length(cell_id)=", length(cell_id), "\n", sep=""), stderr())
	q()
}

seurat_cell_id <- gsub(".*_", "", colnames(x = object))
meta <- rep(0.0, length(seurat_cell_id))
v <- match(seurat_cell_id, cell_id)
meta[seurat_cell_id %in% cell_id] <- values[v[!is.na(v)]]

df <- data.frame(values = meta)
colnames(df) <- gsub("-", "", META_FIELD)
rownames(df) <- colnames(x = object)
object <- AddMetaData(object = object, metadata = df, col.name = colnames(df))

save(object, file=OBJECT)

sessionInfo()

q()

