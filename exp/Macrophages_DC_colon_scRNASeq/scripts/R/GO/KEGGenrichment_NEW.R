
# Perform gene annotation from gene SYMBOL using org.Mm.eg.db database

organism <- "Mm"
orgname <- "mouse"

pval <- 0.01
qval <- 0.05

cat_to_show <- 15

TITLE <- ""

# ask for input
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
	cat("Usage: <genes_list> <out_prefix> [<title>]\n")
	cat("<genes_list>     file with the list of gene SYMBOLS and group, one per line, with header\n")
	cat("                 the group is typically \"upregulated\" or \"downregulated\"\n")
	cat("<out_prefix>     where .tsv table with PW enrichment analysis output is stored\n")
	q() 
}

DATA <- args[1]
OUT_PREFIX <- args[2]
if (length(args) > 2) {
	TITLE <- args[3]
}

orglib <- paste("org.", organism, ".eg.db", sep="")

library(orglib, character.only = TRUE)
library("clusterProfiler")
library("GO.db")
library("limma") # for alias2Symbol()

# See also: https://bioconductor.org/packages/release/bioc/vignettes/ReactomePA/inst/doc/ReactomePA.html

data <- read.table(DATA, header=TRUE)

# extract gene names
data[,1] <- as.character(data[,1])
data[,1] <- alias2SymbolTable(data[,1], species = organism) 
data <- data[!is.na(data[,1]),]

SymbToEnt <- bitr(data[,1], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = orglib, drop=FALSE)

# prepare the data frame
df <- data.frame(Entrez = SymbToEnt$ENTREZID[!is.na(SymbToEnt$ENTREZID)], group = data[,2])

# perform pathway enrichment analysis
pwe <- compareCluster(Entrez~group, data=df, fun="enrichKEGG", organism=orgname, pvalueCutoff = pval, pAdjustMethod = "BH", qvalueCutoff = qval)

if (length(pwe@compareClusterResult$Description) > 0) {

	write.table(as.data.frame(pwe), file = paste(OUT_PREFIX, ".tsv", sep=""), sep="\t", quote=FALSE, row.names=FALSE)

	PLOT <- paste(OUT_PREFIX, ".pdf", sep="")
	height <- max(4,0.5+0.4*min(length(pwe@compareClusterResult$Description), cat_to_show))
	title <- TITLE
	if (TITLE != "") {
		title <- paste(title, " - ", sep=" ")
	} 
	title <- paste(title, "PW enrichment", sep="")
	pdf(PLOT, width=13, height=height)
	print(dotplot(pwe, showCategory=cat_to_show, title=title))
	dev.off()
}


q()


