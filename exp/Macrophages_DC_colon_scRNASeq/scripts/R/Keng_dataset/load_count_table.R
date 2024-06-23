
# Created on 01/16/2020

k <- 20

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
	cat("\nUsage: Rscript load_count_table.R <IN_DIR> <OUT_DIR>\n\n")
	cat("\n<IN_DIR>    where gene-cell tables are stored\n")
	cat("<OUT_DIR>   where objects are saved\n")
	q()
}

IN_DIR <- args[1]
OUT_DIR <- args[2]

TABLE1_FILE <- file.path(IN_DIR, "GSE137927_gf_colon_phagocytes_Dropseq_UMI_raw_counts.txt")
TABLE2_FILE <- file.path(IN_DIR, "GSE137927_spf_colon_phagocytes_Dropseq_UMI_raw_counts.txt")
OBJECT1_FILE <- file.path(OUT_DIR, "object_gf.Robj")
OBJECT2_FILE <- file.path(OUT_DIR, "object_spf.Robj")

LIB_PATH <- file.path(.libPaths()[1], "Seurat_3.0.0")
library("Seurat", lib.loc = LIB_PATH)
library("Matrix")

data1 <- read.table(TABLE1_FILE)
M1 <- Matrix(as.matrix(data1))
object1 <- CreateSeuratObject(counts = M1)
save(object1, file = OBJECT1_FILE)

data2 <- read.table(TABLE2_FILE)
M2 <- Matrix(as.matrix(data2))
object2 <- CreateSeuratObject(counts = M2)
save(object2, file = OBJECT2_FILE)

q()

