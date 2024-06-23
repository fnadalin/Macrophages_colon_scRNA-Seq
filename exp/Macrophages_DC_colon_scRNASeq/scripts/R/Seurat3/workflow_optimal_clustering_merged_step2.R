#!/usr/bin/R 

###### SEURAT V3 #######

THREADS <- 8

SilhouetteSingle <- function(dist, ident, out.prefix = "silhouette") {

	id <- unique(ident)
	if (length(id) > 1) { # FIXED: without this check, it fails if there is only one cluster
		silh <- silhouette(as.numeric(ident), dist)
		summ <- summary(silh)
		save(summ, file = paste(out.prefix, "Robj", sep="."))
		pdf(paste(out.prefix, "pdf", sep="."))
		plot(silh)
		dev.off()
	}

	return()
}


args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
	cat("\nUsage: <obj.Robj> <in_dir> <out_dir> <params>\n")
	cat("\n<obj.Robj>   input Seurat object created at step 1\n")
	cat("<in_dir>     input directory containing the number of PCs and the list of features computed at step 1\n")
	cat("<out_dir>    output directory containing the results of clusters validation\n")
	cat("<params>     file including the values for the parameters, separated by \"=\"\n\n")
	q()
}

LIB_PATH <- file.path(.libPaths()[1], "Seurat_3.0.0")
library("Seurat", lib.loc = LIB_PATH)
library("cluster")
library("parallelDist")

OBJECT <- args[1]
IN_DIR <- args[2]
OUT_DIR <- args[3]
PARAMS <- args[4]

params <- as.matrix(read.table(PARAMS, sep="="))
disps <- as.numeric(unlist(strsplit(params[params[,1]=="disps", 2], split=",")))
nfeats <- as.numeric(unlist(strsplit(params[params[,1]=="nfeats", 2], split=",")))
k <- as.numeric(unlist(strsplit(params[params[,1]=="k", 2], split=",")))
res <- as.numeric(unlist(strsplit(params[params[,1]=="res", 2], split=",")))

object <- eval(parse(text=load(OBJECT)))

cases <- c()
for (ymin in disps)
	cases <- c(cases, paste("mean.var.plot_disp", ymin, sep=""))
for (nfeat in nfeats)
	cases <- c(cases, paste("vst_top", nfeat, sep=""))
drs <- paste("pca", cases, sep="_")

for (ll in c(cases, drs)) 
	dir.create(file.path(OUT_DIR, ll), showWarnings = FALSE, recursive = TRUE)

for (case in cases) {
	print(paste("features =", case))
	features_file <- file.path(IN_DIR, case, "features.txt")
	features <- c(as.matrix(read.table(features_file)))
	data <- t(as.matrix(GetAssayData(object)[features,]))
	dist <- parDist(data, threads = THREADS)
	for (kk in k) {
		for (r in res) {
			cl_ident <- paste("clusters_", case, "_k", kk, "_res", r, sep="")
			out_prefix <- file.path(OUT_DIR, case, cl_ident)
			SilhouetteSingle(dist = dist, ident = c(as.matrix(object@meta.data[cl_ident])), out.prefix = out_prefix) 
		}
	}
}

i <- 1
for (dr in drs) {
	print(paste("dr =", dr))
	case <- cases[i]
	dim_file <- file.path(IN_DIR, case, "num_PCs.txt")
	dim <- c(as.matrix(read.table(dim_file)))
	data <- as.matrix(object@reductions[[dr]]@cell.embeddings[,1:dim])
	dist <- parDist(data, threads = THREADS)
	for (kk in k) {
		for (r in res) {
			cl_ident <- paste("clusters_", dr, "_k", kk, "_res", r, sep="")
			out_prefix <- file.path(OUT_DIR, dr, cl_ident)
			SilhouetteSingle(dist = dist, ident = c(as.matrix(object@meta.data[cl_ident])), out.prefix = out_prefix) 
		}
	}
	i <- i + 1
}


q()

