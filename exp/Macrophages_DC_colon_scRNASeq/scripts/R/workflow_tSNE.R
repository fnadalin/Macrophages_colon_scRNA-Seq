# Global parameters

# seed for tSNE
seed <- 1:10

# perplexity for tSNE
perplexity <- 5*(1:4)

res <- NULL

pt.size <- 1


args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
	cat("\nUsage: Rscript workflow_tSNE.R <object> <out_dir> <mode> [<dim> <seed> <perplexity> <res>]\n\n")
	cat("<object>        .Robj file containing Seurat object \"obj_norm\"\n")
	cat("<out_dir>       where output will be stored\n")
	cat("<mode>          for clustering:\n")
	cat("                \"pc_mouse\": PCs computed on all mouse genes\n")
	cat("                \"hvg\": highly variable mouse genes\n")
	cat("                \"pc_mouse_hvg\": PCs computed on highly variable mouse genes\n")
	cat(paste("<dim>           number of PC to consider (only used in \"pc_notHIV\" or \"pc_hvg\" mode)\n", sep=""))
	cat(paste("<seed>          seeds for the tSNE (default: ", deparse(seed), ")\n", sep=""))
	cat(paste("<perplexity>    tSNE perplexity parameters (default: ", deparse(perplexity), ")\n\n", sep=""))
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
		perplexity <- eval(parse(text=args[6]))
		if (length(args) > 6) {
			res <- eval(parse(text=args[7]))
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
sink(file = paste(OUT_DIR, "/parameters_tsne.txt", sep = ""))
cat(paste("OBJECT=", OBJECT, "\n", sep=""))
cat(paste("OUT_DIR=", OUT_DIR, "\n", sep=""))
cat(paste("mode=", mode, "\n", sep=""))
cat(paste("dim=", dim, "\n", sep=""))
cat(paste("seed=", deparse(seed), "\n", sep=""))
cat(paste("perplexity=", deparse(perplexity), "\n", sep=""))
cat(paste("res=", res, "\n", sep=""))
cat(paste("pt.size=", pt.size, "\n", sep=""))
sink()


object <- LoadObject(OBJECT)

prefix <- mode 
if (dim > 0) {
	prefix <- paste(prefix, dim, sep=".") 
}
reduction <- paste("tsne.", prefix, sep="")

out_dir_tsne <- paste(OUT_DIR, reduction, sep="/")
dir.create(out_dir_tsne, showWarnings = FALSE)

for (s in seed) {
	for (p in perplexity) {
		red <- paste(reduction, ".seed", s, ".perp", p, sep="")	
		cat(paste("Compute ", red, "...", sep=""))
		if (mode == "pc_mouse") {	
			# Compute the tSNE from PCs
			object <- RunTSNE(object, do.fast = TRUE, reduction.use = "pca_mouse", dims.use = 1:dim, seed.use = s, reduction.name = red, perplexity = p)
		} else {
			if (mode == "hvg") {
				# Compute the tSNE from highly variable genes
				hv_genes <- object@var.genes
				object <- RunTSNE(object, do.fast = TRUE, genes.use = hv_genes, seed.use = s, reduction.name = red, perplexity = p)
			} else { # mode == "pc_mouse_hvg"
				# Compute the tSNE from PCs
				object <- RunTSNE(object, do.fast = TRUE, reduction.use = "pca_mouse_hvg", dims.use = 1:dim, seed.use = s, reduction.name = red, perplexity = p)
			}
		}
		cat("done\n")
	}
}

save(object, file = paste(OUT_DIR, "/object.Robj", sep=""))

tSNEplotSamplesMulti(object, out_dir_tsne, seed, perplexity, pt.size = pt.size, reduction.use = reduction)
if (!is.null(res)) {
	tSNEplotClustersMulti(object, out_dir_tsne, seed, perplexity, res = res, prefix = prefix, pt.size = pt.size, reduction.use = reduction)
}



q()


