# Created on 01/27/2020

# Starting from two gene lists, compute the number of ligand-receptor (LR) and receptor-ligand (RL)
# interactions annotated in LRBase.Mmu.eg.db database

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
	cat("\nUsage: <genes1> <genes2> <out>\n")
	cat("\n<genes1>     genes for set 1, one per line\n")
	cat("<genes2>     genes for set 2, one per line\n")
	cat("<out>        output file with info\n")
	q()
}

library("LRBase.Mmu.eg.db")
library("clusterProfiler")
library("limma")

GENES1 <- args[1]
GENES2 <- args[2]
OUT <- args[3]

res <- try(read.table(GENES1))
if (class(res) == 'try-error') {
	genes1 <- c()
} else {
	genes1 <- drop(as.matrix(read.table(GENES1)))
} 

res <- try(read.table(GENES2))
if (class(res) == 'try-error') {
	genes2 <- c()
} else {
	genes2 <- drop(as.matrix(read.table(GENES2)))
} 

genes1_official <- unique(alias2SymbolTable(genes1, species = "Mm"))
genes2_official <- unique(alias2SymbolTable(genes2, species = "Mm"))

genes1_Symb2Ent <- bitr(genes1_official, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db", drop=TRUE)
genes2_Symb2Ent <- bitr(genes2_official, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Mm.eg.db", drop=TRUE)

info_L1 <- select(LRBase.Mmu.eg.db, keys = genes1_Symb2Ent$ENTREZID, columns=c('GENEID_L', 'GENEID_R'), keytype='GENEID_L')
info_R1 <- select(LRBase.Mmu.eg.db, keys = genes1_Symb2Ent$ENTREZID, columns=c('GENEID_L', 'GENEID_R'), keytype='GENEID_R')
info_L2 <- select(LRBase.Mmu.eg.db, keys = genes2_Symb2Ent$ENTREZID, columns=c('GENEID_L', 'GENEID_R'), keytype='GENEID_L')
info_R2 <- select(LRBase.Mmu.eg.db, keys = genes2_Symb2Ent$ENTREZID, columns=c('GENEID_L', 'GENEID_R'), keytype='GENEID_R')

info_LR <- info_L1[info_L1$GENEID_R %in% genes2_Symb2Ent$ENTREZID,]
info_RL <- info_L2[info_L2$GENEID_R %in% genes1_Symb2Ent$ENTREZID,]

sink(file = OUT, append = FALSE)
sink()
sink(file = OUT, append = TRUE)
cat(paste("genes1_file=", GENES1, "\n", sep=""))
cat(paste("genes2_file=", GENES2, "\n", sep=""))
cat(paste("num_genes1_original=", length(genes1), "\n", sep=""))
cat(paste("num_genes2_original=", length(genes2), "\n", sep=""))
cat(paste("num_genes1_official=", length(genes1_official), "\n", sep=""))
cat(paste("num_genes2_official=", length(genes2_official), "\n", sep=""))
cat(paste("num_genes1_entrez=", nrow(genes1_Symb2Ent), "\n", sep=""))
cat(paste("num_genes2_entrez=", nrow(genes2_Symb2Ent), "\n", sep=""))
cat(paste("num_genes1_L=", length(unique(info_L1$GENEID_L)), "\n", sep=""))
cat(paste("num_genes1_R=", length(unique(info_R1$GENEID_R)), "\n", sep=""))
cat(paste("num_genes2_L=", length(unique(info_L2$GENEID_L)), "\n", sep=""))
cat(paste("num_genes2_R=", length(unique(info_R2$GENEID_R)), "\n", sep=""))
cat(paste("LR_pairs=", nrow(info_LR), "\n", sep=""))
cat(paste("RL_pairs=", nrow(info_RL), "\n", sep=""))
sink()


q()

