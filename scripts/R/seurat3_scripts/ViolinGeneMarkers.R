# Created on 2019/06/21

IDENT_SLOT <- ""

# ask for input
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
	cat("\nUsage: Rscript ViolinGeneMarkers.R <obj.Robj> <gene_list> <out.pdf> <height> [<ident_slot>]\n")
	cat("\n<obj.Robj>       input Seurat 3 object\n")
	cat("<genes_list>     comma-separated list of gene names to be plot\n")
	cat("<out.pdf>        output violin plot\n")
	cat("<ident_slot>     identity slot to use\n\n")
	q() 
}

PATH="/data/users/fnadalin/curie/scripts/R/seurat3_scripts"
source(paste(PATH, "/functions.R", sep=""))

OBJECT <- args[1]
GENES_LIST <- args[2]
OUT_PDF <- args[3]
HEIGHT <- as.numeric(args[4])
if (length(args) > 4) 
	IDENT_SLOT <- args[5]

genes_list <- unlist(strsplit(GENES_LIST, split = ","))
object <- eval(parse(text=load(OBJECT)))
if (IDENT_SLOT != "")
	Idents(object) <- IDENT_SLOT

genes_list <- genes_list[genes_list %in% rownames(GetAssayData(object))]

if (length(genes_list) == 0) {
	cat("No gene in the list is found in the matrix\n")
	q()
}

NCOL <- 3

NCOL <- min(length(genes_list),NCOL)
NROW <- ceiling(length(genes_list)/NCOL)
pdf(OUT_PDF, width=max(3.7,3*min(length(genes_list),NCOL)), height=HEIGHT*NROW)
VlnPlot(object = object, features = genes_list, pt.size = 0, ncol = NCOL)
dev.off()

q()

