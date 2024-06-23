
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
	cat("\nUsage: Rscript RunPWenrich.R <deg.table> <filt.table> <db> <title> <out_prefix>\n\n")
	cat("<deg.table>     input .tsv file containing all DEGs\n")
	cat("<filt.table>    output .tsv file containing gene names corresponding to filtered DEGs\n")
	cat("<collection>    GSEA collection: \"h\", \"c2\", \"c5\", \"c7\"\n\n")
	q()
}

deg.table <- args[1]
GSEA.table <- args[2] 
collection <- args[3]
title <- args[4]
out_prefix <- args[5]

PATH="/data/users/fnadalin/curie/scripts/R/annotation/functions.R"
source(PATH)

GSEAPreprocessingLogFC(deg.table, GSEA.table, adjpval.filt = 0.05, revert = FALSE)
RunGSEA(table = GSEA.table, out_prefix, collection, org = "Mm", order = "d", title = title, pval = 0.05, adj_pval = 0.1, nPerm = 1000, minGSSize = 10, maxGSSize = 500)

q()

