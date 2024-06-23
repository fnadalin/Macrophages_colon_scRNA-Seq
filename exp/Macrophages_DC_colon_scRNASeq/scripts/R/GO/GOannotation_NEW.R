
# Perform gene annotation from gene SYMBOL using mm.Hs.eg.db database

organism <- "Mm"

classes <- c("CC", "MF", "BP")
descriptions <- c("Cellular component", "Molecular function", "Biological pathway")
levels <- 2:3

cat_to_show <- 30

# ask for input
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
	cat("\nUsage: <genes_list> <out_dir> [<levels>]\n\n")
	cat("<genes_list>     file with the list of gene SYMBOLS, one per line\n")
	cat("<out_dir>        where .tsv table with GO analysis output is stored\n")
	cat(paste("<levels>         GO levels to consider (default: ", levels, ")\n\n", sep=""))
	q() 
}

DATA <- args[1]
OUT_DIR <- args[2]
if (length(args) > 2) {
	levels <- eval(parse(text=args[3]))
}

orglib <- paste("org.", organism, ".eg.db", sep="")

library(orglib, character.only = TRUE)
library("clusterProfiler")
library("GO.db")
library("limma") # for alias2Symbol()

# See also: https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.html#selecting-a-biomart-database-and-dataset
# See also: http://www.pangloss.com/wiki/R/BioMart
# See also: https://www.biostars.org/p/70821/

data <- read.table(DATA)

# extract gene names
names <- c(as.matrix(data))
names <- alias2SymbolTable(names, species = organism)
names <- names[!is.na(names)]

setwd(OUT_DIR)

SymbToEnt <- bitr(names, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = orglib)
# write.csv(x=SymbToEnt, file="SymbToEntrez.csv")

i <- 1
for (class in classes) {
	for (level in levels) {
		# perform gene ontology annotation
		# specify: type of ontology ("MF", "BP", "CC")
		#          level
		ggo <- groupGO(gene = SymbToEnt$ENTREZID, OrgDb = orglib, ont = class, level = level, readable = TRUE)
		file.name <- paste("GO_", class, level, ".tsv", sep="")
		write.table(x=ggo, file=file.name, sep="\t", quote=FALSE)
		
		# Barplots of GO
		plot.name <- paste("barplot_groupGO_", class, "_", level, ".pdf", sep="")
		height <- max(4,0.5+0.4*min(length(ggo@result$Description), cat_to_show))
		TITLE <- paste("GO annotation (", class, ", level=", level, ")", sep="")
		pdf(plot.name, width=9, height=height)
		print(barplot(ggo, drop=TRUE, showCategory=cat_to_show, title=TITLE))
		dev.off()
	}
	i <- i + 1
}

q()


