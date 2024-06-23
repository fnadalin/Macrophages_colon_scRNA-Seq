
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
	cat("\nUsage: Rscript RunPWenrich.R <deg.table> <filt.table> <db> <title> <out_prefix>\n\n")
	cat("<deg.table>     input .tsv file containing all DEGs\n")
	cat("<filt.table>    output .tsv file containing gene names corresponding to filtered DEGs\n")
	cat("<db>            PW database (\"reactome\" or \"kegg\")\n\n")
	q()
}

deg.table <- args[1]
filt.table <- args[2] 
db <- args[3]
title <- args[4]
out_prefix <- args[5]

PATH="/data/users/fnadalin/curie/scripts/R/annotation/functions.R"
source(PATH)

GOenrichPreprocessing(deg.table, filt.table, min.perc = 0.1, logFC.filt = 0.5, adjpval.filt = 0.05)
PWenrich(filt.table, out_prefix, db = db, org = "Mm", orgname = "mouse", pval = 0.5, qval = 1, cat.to.show = 15, title = title) 

q()

