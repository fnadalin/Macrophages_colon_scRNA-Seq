
# Created on 11/28/2019

INFTY <- 30000

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
	cat("\nUsage: Rscript filterDEGs.R <table> <filt_table> <log2FC> <pval> <perc>\n")
	cat("\n<table>       input table containing the DEGs\n")
	cat("<table_filt>  output table containing the filtered DEGs\n")
	cat("<log2FC>      minimum abs log2FC\n")
	cat("<pval>        maximum adj. p-val\n")
	cat("<perc>        minimum percentage of detection (ranging between 0 and 1)\n\n")
	q()
}


TABLE <- args[1]
TABLE_FILT <- args[2]
LOG2FC <- as.numeric(args[3])
PVAL <- as.numeric(args[4])
PERC <- as.numeric(args[5])
library("funr")

WORKING_DIR <- getwd()
SCRIPT_PATH <- dirname(sys.script())
SCRIPT_NAME <- basename(sys.script())
setwd(SCRIPT_PATH)
FUNCTIONS <- "functions.R"
source(FUNCTIONS)
setwd(WORKING_DIR)

deg <- read.table(TABLE, sep="\t", header=TRUE)
deg.filt <- FilterDEGtable(deg, logFC.filt = LOG2FC, adjpval.filt = PVAL, min.pct = PERC, num = INFTY)

write.table(deg.filt, file = TABLE_FILT, sep="\t", quote=FALSE, row.names=FALSE)

q()


