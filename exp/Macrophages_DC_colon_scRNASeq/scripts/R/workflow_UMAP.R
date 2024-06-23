# Global parameters

# seed for UMAP
seed <- 1:10

# neighbors for UMAP
neighbors <- 30

# min dist for UMAP
min_dist <- 0.3

res <- NULL

pt.size <- 1


args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
	cat("\nUsage: Rscript workflow_UMAP.R <object> <out_dir> <mode> [<dim> <seed> <neighbors> <min_dist> <res>]\n\n")
	cat("<object>        .Robj file containing Seurat object \"obj_norm\"\n")
	cat("<out_dir>       where output will be stored\n")
	cat("<mode>          for clustering:\n")
	cat("                \"pc_mouse\": PCs computed on all mouse genes\n")
	cat("                \"hvg\": highly variable mouse genes\n")
	cat("                \"pc_mouse_hvg\": PCs computed on highly variable mouse genes\n")
	cat(paste("<dim>           number of PC to consider (only used in \"pc_notHIV\" or \"pc_hvg\" mode)\n", sep=""))
	cat(paste("<seed>          seeds for the tSNE (default: ", deparse(seed), ")\n", sep=""))
	cat(paste("<neighbors>     UMAP neighbor parameter values (default: ", deparse(neighbors), ")\n\n", sep=""))
	cat(paste("<min_dist>      UMAP min dist parameter values (default: ", deparse(min_dist), ")\n\n", sep=""))
	cat(paste("<res>           resolution parameter for clustering (OPTIONAL)\n", sep=""))
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
	seed <- eval(parse(text=args[5]))
	if (length(args) > 5) {
		neighbors <- eval(parse(text=args[6]))
		if (length(args) > 6) {
			min_dist <- eval(parse(text=args[7]))
			if (length(args) > 7) {
				res <- eval(parse(text=args[8]))
			}
		}
	}
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
sink(file = paste(OUT_DIR, "/parameters_umap.txt", sep = ""))
cat(paste("object=", OBJECT, "\n", sep=""))
cat(paste("out.dir=", OUT_DIR, "\n", sep=""))
cat(paste("mode=", mode, "\n", sep=""))
cat(paste("dim=", dim, "\n", sep=""))
cat(paste("seed=", deparse(seed), "\n", sep=""))
cat(paste("neighbors=", deparse(neighbors), "\n", sep=""))
cat(paste("min_dist=", deparse(min_dist), "\n", sep=""))
cat(paste("res=", res, "\n", sep=""))
cat(paste("pt.size=", pt.size, "\n", sep=""))
sink()


object <- LoadObject(OBJECT)

prefix <- mode 
if (dim > 0) {
	prefix <- paste(prefix, dim, sep=".") 
}
reduction <- paste("umap.", prefix, sep="")

out_dir_tsne <- paste(OUT_DIR, reduction, sep="/")
dir.create(out_dir_tsne, showWarnings = FALSE)

for (s in seed) {
	for (n in neighbors) {
		for (d in min_dist) {
			red <- paste(reduction, ".seed", s, ".neigh", n, ".dist", d, sep="")	
			cat(paste("Compute ", red, "...", sep=""))
			if (mode == "pc_mouse") {	
				# Compute the tSNE from PCs
				hv_genes <- object@var.genes
				object <- RunUMAP(object, do.fast = TRUE, reduction.use = "pca_mouse", dims.use = 1:dim, seed.use = s, reduction.name = red, n_neighbors = n, min_dist = d)
			} else {
				if (mode == "hvg") {
					# Compute the tSNE from highly variable genes
					hv_genes <- object@var.genes
					object <- RunUMAP(object, do.fast = TRUE, genes.use = hv_genes, seed.use = s, reduction.name = red, n_neighbors = n, min_dist = d)
				} else { # mode == "pc_mouse_hvg"
					# Compute the tSNE from PCs
					object <- RunUMAP(object, do.fast = TRUE, reduction.use = "pca_mouse_hvg", dims.use = 1:dim, seed.use = s, reduction.name = red, n_neighbors = n, min_dist = d)
				}
			}
			cat("done\n")
		}
	}
}

save(object, file = paste(OUT_DIR, "/object.Robj", sep=""))

UMAPplotSamplesMulti(object, out_dir_tsne, seed, neighbors = neighbors, min_dist = min_dist, pt.size = pt.size, reduction.use = reduction)
if (!is.null(res)) {
	UMAPplotClustersMulti(object, out_dir_tsne, seed, res = res, prefix = prefix, neighbors = neighbors, min_dist = min_dist, pt.size = pt.size, reduction.use = reduction)
}



q()


