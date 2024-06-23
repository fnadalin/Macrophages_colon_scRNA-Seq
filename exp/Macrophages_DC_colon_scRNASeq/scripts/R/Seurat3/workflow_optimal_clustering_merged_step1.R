#!/usr/bin/R 

###### SEURAT V3 #######

# 1. Select the top significant PCs
# 2. Run the clustering with different parameters (number of neighbors and resolution)


LIB_PATH <- file.path(.libPaths()[1], "Seurat_3.0.0")
library("Seurat", lib.loc = LIB_PATH)
library("future")


ClustDr <- function(object, dr = "pca", pcs = 50, pval = 1e-05, non.rand.sd.frac = 0.25, k, res, file) {

	# compute the number of top significant PCs 
	# that explain 25% more variance compared to random (i.e. last 10 PCs)
	# this check is meant to exclude artifacts due to the svd approximation
	# that is used to compute the PCs
	min_stdev <- (1+non.rand.sd.frac)*mean(object@reductions[[dr]]@stdev[(pcs-10):pcs])
	nPCs <- 0
	for (p in object@reductions[[dr]]@jackstraw@overall.p.values[,2]) {
		if (p > pval || object@reductions[[dr]]@stdev[nPCs+1] < min_stdev) 
			break
		nPCs <- nPCs + 1
	}
	# take the highest PC that explains > 0.05 variance with respect to the next
#	max_sd_PC <- max(which(sd[1:(pcs-1)]-sd[2:pcs] > delta_sd))
	write(nPCs, file = file)

	for (kk in k) {
		object <- FindNeighbors(object, k.param = kk, force.recalc = TRUE, reduction = dr, dims = 1:nPCs)
		for (r in res) {
			cl_ident <- paste("clusters_", dr, "_k", kk, "_res", r, sep="")
			object <- FindClusters(object, resolution = r)
			names(object@meta.data)[names(object@meta.data) == "seurat_clusters"] <- cl_ident
			object@active.ident <- as.factor(object$orig.ident)
		}
	}

	return(object)
}


ClustGenes <- function(object, genes, genes_list_label = "hvg", k, res) {

	for (kk in k) {
		object <- FindNeighbors(object, k.param = kk, force.recalc = TRUE, features = genes, dims = NULL)
		for (r in res) {
			cl_ident <- paste("clusters_", genes_list_label, "_k", kk, "_res", r, sep="")
			object <- FindClusters(object, resolution = r)
			names(object@meta.data)[names(object@meta.data) == "seurat_clusters"] <- cl_ident
			object@active.ident <- as.factor(object$orig.ident)
		}
	}

	return(object)
}




args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
	cat("\nUsage: <obj.Robj> <out_dir> <params>\n")
	cat("\n<obj.Robj>   input Seurat object created at step 0\n")
	cat("<out_dir>    output directory containing the list of genes (from previous step) and the results of PCs selection\n")
	cat("<params>     file including the values for the parameters, separated by \"=\"\n\n")
	q()
}

OBJECT <- args[1]
OUT_DIR <- args[2]
PARAMS <- args[3]


params <- as.matrix(read.table(PARAMS, sep="="))
cores <- as.numeric(params[params[,1]=="cores", 2])
disps <- as.numeric(unlist(strsplit(params[params[,1]=="disps", 2], split=",")))
nfeats <- as.numeric(unlist(strsplit(params[params[,1]=="nfeats", 2], split=",")))
k <- as.numeric(unlist(strsplit(params[params[,1]=="k", 2], split=",")))
res <- as.numeric(unlist(strsplit(params[params[,1]=="res", 2], split=",")))


plan("multiprocess", workers = cores)

object <- eval(parse(text=load(OBJECT)))

for (feature_method in c("mean.var.plot", "vst")) {
	if (feature_method == "mean.var.plot") {
		for (ymin in disps) {
			case <- paste("mean.var.plot_disp", ymin, sep="")
			dr <- paste("pca", case, sep="_")
			out_dir <- file.path(OUT_DIR, case)
			object <- ClustDr(object, dr = dr, k = k, res = res, file = file.path(out_dir, "num_PCs.txt"))
			genes_file <- file.path(OUT_DIR, case, "features.txt")
			object <- ClustGenes(object, genes = c(as.matrix(read.table(genes_file))), genes_list_label = case, k = k, res = res)		
		}
	} else { 
		for (nfeat in nfeats) {
			case <- paste("vst_top", nfeat, sep="")
			dr <- paste("pca", case, sep="_")
			out_dir <- file.path(OUT_DIR, case)
			object <- ClustDr(object, dr = dr, k = k, res = res, file = file.path(out_dir, "num_PCs.txt"))
			genes_file <- file.path(OUT_DIR, case, "features.txt")
			object <- ClustGenes(object, genes = c(as.matrix(read.table(genes_file))), genes_list_label = case, k = k, res = res)	
		}
	}
}

save(object, file=OBJECT)



q()

