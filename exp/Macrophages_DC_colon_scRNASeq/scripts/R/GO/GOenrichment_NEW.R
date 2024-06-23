
# Perform gene annotation from gene SYMBOL using org.Hs.eg.db database
# See also: http://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html

organism <- "Mm"

pval <- 0.01
qval <- 0.05
LEVEL <- 4 # granularity level for GO terms filtering
SIMILARITY <- 0.7 # similarity cut-off between GO terms

cat_to_show <- 15

TITLE <- ""

# ask for input
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
	cat("Usage: <genes_list> <out_prefix> [<title>]\n")
	cat("<genes_list>     file with the list of gene SYMBOLS and value, one per line, with header\n")
	cat("                 the value is the name of the group, typically \"upregulated\" or \"downregulated\"\n")
	cat("<out_prefix>     where .tsv table with GO enrichment analysis output is stored\n")
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

data <- read.table(DATA, header=TRUE, sep="\t")

# extract gene names
data[,1] <- as.character(data[,1])
data[,1] <- alias2SymbolTable(data[,1], species = organism) 
data <- data[!is.na(data[,1]),]

SymbToEnt <- bitr(data[,1], fromType = "SYMBOL", toType = "ENTREZID", OrgDb = orglib, drop=FALSE)

# prepare the data frame
df <- data.frame(Entrez = SymbToEnt$ENTREZID[!is.na(SymbToEnt$ENTREZID)], group = data[,2])

# perform gene ontology enrichment analysis
# specify: type of ontology ("MF", "BP", "CC")
for (type in c("MF", "BP", "CC")) {
	# perform GO enrichment analysis
	ego <- compareCluster(Entrez~group, data=df, fun="enrichGO", OrgDb = orglib, ont = type, pAdjustMethod = "BH", pvalueCutoff  = pval, qvalueCutoff = qval, readable = TRUE)

	if (length(ego@compareClusterResult$Description)  > 0) {
	
		# set the minimum granularity level
		#	ego <- gofilter(ego, level = LEVEL)
		# remove redundant GO terms and by keeping a representative terms
		# the representative term is the one with lowest p-value
		ego <- simplify(ego, cutoff = SIMILARITY, by = "p.adjust", select_fun = min)

		write.table(as.data.frame(ego), file = paste(OUT_PREFIX, "_", type, ".tsv", sep=""), sep="\t", quote=FALSE, row.names=FALSE)
		
		PLOT <- paste(OUT_PREFIX, "_", type, ".pdf", sep="")
		height <- max(4,0.5+0.4*min(length(ego@compareClusterResult$Description), cat_to_show))
		title <- TITLE
		if (TITLE != "") {
			title <- paste(title, " - ", sep=" ")
		} 
		title <- paste(title, "GO enrichment", sep="")
		pdf(PLOT, width=12, height=height)
		print(dotplot(ego, showCategory=cat_to_show, title=title))
		dev.off()
	}
}


q()


