#!/usr/bin/R 

###### SEURAT V3 #######


args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
	cat("\nUsage: <obj.Robj> <in_dir> <out_dir> <out_table> <params>\n\n")
	cat("\n<obj.Robj>   input Seurat object created at step 2\n")
	cat("<in_dir>     input directory containing the results of silhouette scoring\n")
	cat("<out_dir>    output directory containing the average expression of gene markers\n")
	cat("<out_table>  output .tsv file where the summary of the clustering validation is saved\n")
	cat("<params>     file including the values for the parameters, separated by \"=\"\n\n")
	q()
}

LIB_PATH <- file.path(.libPaths()[1], "Seurat_3.0.0")
library("Seurat", lib.loc = LIB_PATH)

OBJECT <- args[1]
IN_DIR <- args[2]
OUT_DIR <- args[3]
OUT_TABLE <- args[4]
PARAMS <- args[5]

params <- as.matrix(read.table(PARAMS, sep="="))
disps <- as.numeric(unlist(strsplit(params[params[,1]=="disps", 2], split=",")))
nfeats <- as.numeric(unlist(strsplit(params[params[,1]=="nfeats", 2], split=",")))
k <- as.numeric(unlist(strsplit(params[params[,1]=="k", 2], split=",")))
res <- as.numeric(unlist(strsplit(params[params[,1]=="res", 2], split=",")))
gene_markers <- unlist(strsplit(params[params[,1]=="gene_markers", 2], split=","))
conditions <- unlist(strsplit(params[params[,1]=="conditions", 2], split=","))

num_genes <- length(gene_markers)

object <- eval(parse(text=load(OBJECT)))

cases <- c()
for (ymin in disps) 
	cases <- c(cases, paste("mean.var.plot_disp", ymin, sep=""))
for (nfeat in nfeats)
	cases <- c(cases, paste("vst_top", nfeat, sep=""))
drs <- paste("pca", cases, sep="_")

# write header
header <- c("features_or_dr", "k")
header <- c(header, "cluster_res", "num_clusters", "BM_cluster_ID", "BM_cluster_size", "BM_distal_frac", "BM_proximal_frac", "avg_silhouette", "BM_avg_silhouette")
for (g in gene_markers)
	header <- c(header, paste("BM_cluster_avg", g, "expr", sep="_"))
write(paste(header, collapse = "	"), file = OUT_TABLE)

for (ll in c(cases,drs)) {
	print(paste("features/dr =", ll))
	dir.create(file.path(OUT_DIR, ll), showWarnings = FALSE, recursive = TRUE)
	for (kk in k) {
		print(paste("k =", kk))
		for (r in res) {
			ID <- paste("clusters_", ll, "_k", kk, "_res", r, sep="")
			silh_file <- file.path(IN_DIR, ll, paste(ID, "Robj", sep="."))
			if (file.exists(silh_file)) {
				object <- SetIdent(object = object, value = ID)
				new_ident <- Idents(object = object)
				cl <- unique(new_ident)
				num_cl <- length(cl)
				expr <- AverageExpression(object = object, features = gene_markers)$RNA
				# print gene markers expression
				expr_file <- file.path(OUT_DIR, ll, paste(ID, "txt", sep="."))
				write.table(expr, file = expr_file, sep = "\t", quote = FALSE, col.names = as.character(colnames(expr)))
				# compute the cluster size 
				cl_size <- unlist(lapply(cl, function(x) sum(new_ident == x)))
				# compute the fraction of cells in each cluster for the 2 conditions
				cells_cond <- object$orig.ident %in% conditions
				frac_in_cond <- unlist(lapply(cl, function(x) sum((new_ident == x) & cells_cond) / sum(cells_cond)))
				frac_not_in_cond <- unlist(lapply(cl, function(x) sum((new_ident == x) & !cells_cond) / sum(!cells_cond)))
				# print cluster composition
				M <- t(matrix(c(cl_size,frac_in_cond,frac_not_in_cond), nrow=num_cl, ncol=3))
				colnames(M) <- cl
				rownames(M) <- c("cl_size", "frac_in_cond", "frac_not_in_cond")
				comp_file <- file.path(OUT_DIR, ll, paste(ID, "_comp.txt", sep=""))
				write.table(M, file = comp_file, sep = "\t", quote = FALSE, col.names = as.character(colnames(M)))
				# detect the biggest cluster that contains > 50% of cells in the condition(s)
				idx <- order(cl_size, decreasing = TRUE)
				cl_ord <- cl[idx] 
				sel_cl_name <- cl_ord[1]
				if (sum(frac_in_cond > frac_not_in_cond) > 0) {
					cl_ord_sel <- cl_ord[frac_in_cond[idx] > frac_not_in_cond[idx]]
					sel_cl_name <- cl_ord_sel[1]
				} 
				num_sel_cl <- length(WhichCells(object, idents = sel_cl_name))
				# parse silhouette summary
				silh <- eval(parse(text = load(silh_file)))
				# print info to table
				values <- c(ll, kk, r, num_cl, as.character(sel_cl_name), num_sel_cl, frac_in_cond[cl == sel_cl_name], frac_not_in_cond[cl == sel_cl_name], silh$avg.width, silh$clus.avg.width[sel_cl_name], expr[,sel_cl_name])
				write(paste(values, collapse = "	"), file = OUT_TABLE, append = TRUE)
			}
		}
	}
}



q()

