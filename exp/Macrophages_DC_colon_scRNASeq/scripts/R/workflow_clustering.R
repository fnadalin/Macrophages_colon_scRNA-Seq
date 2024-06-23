# Global parameters

# resolution parameter for clustering
res <- 0.1*(1:8)
dim <- 0

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
	cat("\nUsage: Rscript workflow_clustering.R <object> <out_dir> <mode> [<dim> <res>]\n\n")
	cat("<object>        .Robj file containing Seurat object \"obj_norm\"\n")
	cat("<out_dir>       where output will be stored\n")
	cat("<mode>          for clustering:\n")
	cat("                \"pc_mouse\": PCs computed on mouse genes\n")
	cat("                \"hvg\": highly variable mouse genes\n")
	cat("                \"pc_mouse_hvg\": PCs computed on highly variable mouse genes\n")
	cat(paste("<dim>           number of PC to consider (only used in \"pc_mouse\" or \"pc_mouse_hvg\" mode)\n", sep=""))
	cat(paste("<res>           resolution parameter for clustering (default: ", deparse(res), ")\n", sep=""))
	q()
}

OBJECT <- args[1]
OUT_DIR <- args[2]
mode <- args[3]

if (!(mode == "pc_mouse" || mode == "hvg" || mode == "pc_mouse_hvg")) {
	cat("<mode> must be either \"pc_mouse\" or \"hvg\" or \"pc_mouse_hvg\"\n")
	q()
}

if (mode == "pc_mouse" || mode == "pc_mouse_hvg") {
	if (length(args) < 4) {
		cat(paste("<dim> parameter is required for mode ", mode, "\n", sep=""))
		q()
	} else {
		dim <- as.numeric(args[4])
	}
}
if (length(args) > 4) {
	res <- eval(parse(text=args[5]))
}


library("funr")

WORKING_DIR <- getwd()
SCRIPT_PATH <- dirname(sys.script())
SCRIPT_NAME <- basename(sys.script())
setwd(SCRIPT_PATH)
FUNCTIONS <- "functions.R"
source(FUNCTIONS)
setwd(WORKING_DIR)


# print the parameters to file
sink(file = paste(OUT_DIR, "/parameters_cl_", mode, ".txt", sep = ""))
cat(paste("OBJECT=", OBJECT, "\n", sep=""))
cat(paste("OUT_DIR=", OUT_DIR, "\n", sep=""))
cat(paste("mode=", mode, "\n", sep=""))
cat(paste("dim=", dim, "\n", sep=""))
cat(paste("res=", deparse(res), "\n", sep=""))
sink()


object <- LoadObject(OBJECT)

# for columns in object@meta.data identifying clusters
prefix <- mode 
if (dim > 0) {
	prefix <- paste(prefix, dim, sep=".") 
}

# Assign cells to clusters and save the identities to id = "res.$r" where r is the resolution parameter (ranging from 0.1 to 0.8 by default)
if (mode == "pc_mouse") {
	# Cluster the cells using a graph algorithm from PCs
	object <- CellsClusters(object, out.dir = OUT_DIR, dim = dim, res = res, prefix = prefix, reduction.type = "pca_mouse", suffix = mode)
	# data <- as.matrix(object@dr$pca@cell.embeddings[1:dim]))
	data <- as.matrix(eval(parse(text = paste("object@dr$pca_mouse@cell.embeddings[,1:", dim, "]", sep=""))))
} else {
	if (mode == "hvg") {
		# Cluster the cells using a graph algorithm from highly variable genes
		hv_genes <- object@var.genes
		obj_norm_hvg <- CellsClusters(object, out.dir = OUT_DIR, dim = dim, res = res, prefix = prefix, genes.use = hv_genes, suffix = mode)
		data <- t(as.matrix(object@var.genes))
	} else { # "pc_mouse_hvg"
		# Cluster the cells using a graph algorithm from PCs
		object <- CellsClusters(object, out.dir = OUT_DIR, dim = dim, res = res, prefix = prefix, reduction.type = "pca_mouse_hvg", suffix = mode)
		# data <- as.matrix(object@dr$pca_hvg@cell.embeddings[,1:dim])
		data <- as.matrix(eval(parse(text = paste("object@dr$pca_mouse_hvg@cell.embeddings[,1:", dim, "]", sep=""))))
	}
}


save(object, file = paste(OUT_DIR, "/object.Robj", sep=""))


# Validate clusters with Silhouette method
out_dir_sil <- paste(OUT_DIR, "/", prefix, sep="")
dir.create(out_dir_sil, showWarnings = FALSE)

out_prefix <- paste(out_dir_sil, "/cl_", prefix, sep = "")
Silhouette(object, data, out.prefix = out_prefix, res, prefix = prefix)




q()


