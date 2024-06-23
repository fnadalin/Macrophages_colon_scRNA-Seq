
# Created on 01/17/2020

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
	cat("\nUsage: Rscript reformat_DEG_table.R <in.tsv> <out.tsv>\n\n")
	q()
}

IN_TSV <- args[1]
OUT_TSV <- args[2]

data <- read.table(IN_TSV, sep = "\t", header=TRUE)

colnames(data)[1] <- "geneID"
data$avg_logFC <- data$avg_logFC/log(2)
colnames(data)[colnames(data) == "avg_logFC"] <- "avg_log2FC"

write.table(data, file = OUT_TSV, quote = FALSE, sep = "\t", row.names = FALSE)

q()
