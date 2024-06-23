
# Perform GSEA from gene SYMBOL using org.Hs.eg.db database

# See also: http://software.broadinstitute.org/gsea/msigdb/collections.jsp#H
MSIGDB_PREFIX <- "~/curie/exp/Macrophages_DC_colon_scRNASeq/data/MSigDb/"
MSIGDB_SUFFIX <- ".all.v6.2.entrez.gmt"

pval <- 0.05
adj_pval <- 0.1
nPerm <- 1000
minGSSize <- 10
maxGSSize <- 500

order <- "a"

TITLE <- ""


# ask for input
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
	cat("\nUsage: Rscript GSEA.R <genes_list> <out_prefix> <c> [<o> <title>]\n")
	cat("<genes_list>     file with the list of gene SYMBOLS and value (t statistic or fold change) [one per line, with header]\n")
	cat("<out_prefix>     where GSEA output is stored\n")
	cat("<c>              MSigDb collection: \"h\", \"c2\", \"c5\", \"c7\"\n")
	cat(paste("<o>              order score in ascending (\"a\") or descending (\"d\") order [OPTIONAL - default: ", order, "]\n\n", sep=""))
	q() 
}

DATA <- args[1]
OUT_PREFIX <- args[2]
COLLECTION <- args[3]
if (length(args) > 3) {
	order <- args[4]
	if (length(args) > 4) {
		TITLE <- args[5]
	}
}

library("org.Mm.eg.db")
library("org.Hs.eg.db")
library("clusterProfiler")
library("GSEABase")
library("limma") # for alias2Symbol()
library("biomaRt")

mart1 = useMart("ensembl", dataset="mmusculus_gene_ensembl")
mart2 = useMart("ensembl", dataset="hsapiens_gene_ensembl") 

# See also: https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.html#selecting-a-biomart-database-and-dataset
# See also: http://www.pangloss.com/wiki/R/BioMartENTREZID
# See also: https://www.biostars.org/p/70821/

data <- read.table(DATA, header = TRUE)

data[,1] <- alias2SymbolTable(data[,1], species = "Mm") # to convert gene aliases to official gene names
data <- data[!is.na(data[,1]),]
SymbToEns_mouse <- bitr(data[,1], fromType = "SYMBOL", toType = "ENSEMBL", OrgDb = org.Mm.eg.db, drop = FALSE)

# mouse -> human
Ens_mouse_human <- getLDS (attributes=c("ensembl_gene_id"),
						   filters="ensembl_gene_id", values=SymbToEns_mouse$ENSEMBL, mart=mart1,
						   attributesL=c("ensembl_gene_id"), martL=mart2)
Ens_human <- Ens_mouse_human$Gene.stable.ID.1

EnsToEnt_human <- bitr(Ens_human, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop = FALSE)
Ent_human <- EnsToEnt_human$ENTREZID


########### Map the indexes #############

# v[i] = index of EnsToEnt_human$ENSEMBL[i] in Ens_mouse_human
v <- match(EnsToEnt_human$ENSEMBL, Ens_mouse_human$Gene.stable.ID.1)
# w[j] = index of Ens_mouse_human$Gene.stable.ID[j] in SymbToEns_mouse
w <- match(Ens_mouse_human$Gene.stable.ID, SymbToEns_mouse$ENSEMBL)
# x[h] = index of SymbToEns_mouse$SYMBOL[h] in data
x <- match(SymbToEns_mouse$SYMBOL, data[,1])
# reorder data based on Ent_human ordering
data_ordered <- data[x[w[v]],2]
names(data_ordered) <- EnsToEnt_human$ENTREZ
## y[k] = index of unique elements in EnsToEnt_human$ENTREZ
y <- match(unique(EnsToEnt_human$ENTREZ),  EnsToEnt_human$ENTREZ)
## select a unique value for each entrez identifier
data_ordered <- data_ordered[y]


####### Prepare the gene list ###########

## feature 1: numeric vector
geneList <- data_ordered[!is.na(data_ordered) && !is.na(names(data_ordered))]
## feature 3: decreasing order
if (order == "a") {
	geneList <- 0-geneList
}
geneList <- sort(geneList, decreasing = TRUE)

geneList <- geneList[!is.na(names(geneList))]


############### Run GSEA ################

# the collection contains ENTREZ gene IDs
msigdb <- paste(MSIGDB_PREFIX, COLLECTION, MSIGDB_SUFFIX, sep="")
collection <- read.gmt(msigdb)

gsea <- GSEA(geneList = geneList, exponent = 1, nPerm = nPerm, minGSSize = minGSSize, maxGSSize = maxGSSize, TERM2GENE = collection, pvalueCutoff = adj_pval, pAdjustMethod = "BH", verbose = FALSE)
# gsea_res <- gsea@result[gsea@result$pvalue < pval,]

write.table(gsea@result, file = paste(OUT_PREFIX, ".tsv", sep=""), sep="\t", quote=FALSE, row.names=FALSE)

sign_gene_sets <- which(gsea@result$p.adjust < adj_pval)
if (length(sign_gene_sets) == 0)
{
	q()
}


######## Plot the enriched sets #########

cat_to_show <- min(nrow(gsea@result), 15)
height <- max(4,0.5+0.4*cat_to_show)
title <- TITLE
if (TITLE != "") {
	title <- paste(title, " - ", sep=" ")
}
title <- paste(title, "GSEA (collection \"", COLLECTION, "\")", sep="")

PLOT <- paste(OUT_PREFIX, ".pdf", sep="")
pdf(PLOT, width=12, height=height)
dotplot(gsea, showCategory=cat_to_show, title=title)
dev.off()

if (cat_to_show > 0) {
	for (i in 1:cat_to_show) {
		PLOT <- paste(OUT_PREFIX, "_" , gsea@result$ID[i], ".pdf", sep="")
		pdf(PLOT, width=7, height=7)
		print(gseaplot(gsea, geneSetID = gsea@result$ID[i], title=gsea@result$Description[i]))
		dev.off()
	}
}


q()


