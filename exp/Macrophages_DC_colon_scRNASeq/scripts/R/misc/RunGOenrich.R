
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
	cat("\nUsage: Rscript RunGOenrich.R <deg.table> <filt.table> <title> <out_prefix>\n\n")
	cat("<deg.table>     input .tsv file containing all DEGs\n")
	cat("<filt.table>    output .tsv file containing gene names corresponding to filtered DEGs\n\n")
	q()
}

deg.table <- args[1]
filt.table <- args[2] 
title <- args[3]
out_prefix <- args[4]

PATH="/data/users/fnadalin/curie/scripts/R/annotation/functions.R"
source(PATH)

GOenrichPreprocessing(deg.table, filt.table, min.perc = 0.1, logFC.filt = 0.5, adjpval.filt = 0.05)
GOenrich(filt.table, out_prefix, org = "Mm", pval = 0.01, qval = 0.05, level = 4, similarity = 0.7, cat.to.show = 15, title = title)

q()

