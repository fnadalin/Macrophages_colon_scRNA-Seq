# NB: THE FIRST LINE WAS SKIPPED IN THE PREVIOUS VERSION!!!!

# ask for input
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
	cat("\nUsage: Rscript geneSynonimToSymbol.R <genes_list> <out_list> <org>\n")
	cat("<genes_list>     input file with the list of gene SYMBOLS, one per line\n")
	cat("<out_list>       output file with the list of representative SYMBOLS\n")
	cat("<org>            organism [\"Hs\" or \"Mm\"]\n")
	q() 
}

IN_GENES <- args[1]
OUT_GENES <- args[2]
ORG <- args[3]

library("limma") # for alias2Symbol()

genes <- read.table(IN_GENES)
genes <- drop(as.matrix(genes))
n <- length(genes)

genes_official <- alias2SymbolTable(genes, species=ORG) # to convert gene aliases to official gene names
genes_official <- genes_official[!is.na(genes_official)]
m <- length(genes_official)

cat(paste(n, "genes read,", m , "gene official names\n"))
write.table(genes_official, file=OUT_GENES, quote=FALSE, row.names=FALSE, col.names=FALSE)

q()

