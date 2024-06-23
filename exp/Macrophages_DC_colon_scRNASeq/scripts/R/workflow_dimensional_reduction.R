# Global parameters

# number of PCs to compute
num_pc <- 30
jack <- (0 > 1)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
	cat("\nUsage: Rscript workflow_dimensional_reduction.R <object> <out_dir> <mode> [<num_pc> <j>]\n\n")
	cat("<object>        .Robj file containing Seurat object \"obj_norm\"\n")
	cat("<out_dir>       where the new object is stored\n")
	cat("<mode>          for dimensional reduction:\n")
	cat("                \"mouse\": PCA on mouse genes only\n")
	cat("                \"mouse_hvg\": PCA on highly variable mouse genes\n")	
	cat(paste("<num_pc>        number of PCs to compute (default: ", num_pc, "\n", sep=""))
	cat("<j>             generate Jackstraw plot (computationally expensive - default: FALSE)\n\n")
	q()
}

OBJECT <- args[1]
OUT_DIR <- args[2]
mode <- args[3]
if (length(args) > 3) {
	num_pc <- as.numeric(args[4])
	if (length(args) > 4) {
		jack <- (args[5] == "j")
	}
}


if (!(mode == "mouse" || mode == "mouse_hvg")) {
	cat("<mode> must be either \"mouse\" or \"mouse_hvg\"\n")
	q()
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
sink(file = paste(OUT_DIR, "/parameters_dr.txt", sep = ""))
cat(paste("OBJECT=", OBJECT, "\n", sep=""))
cat(paste("OUT_DIR=", OUT_DIR, "\n", sep=""))
cat(paste("mode=", mode, "\n", sep=""))
cat(paste("num_pc=", num_pc, "\n", sep=""))
cat(paste("jack=", jack, "\n", sep=""))
sink()


cat("Load the Seurat object\n")

object <- LoadObject(OBJECT)

if (mode == "mouse") {
	# Compute the PCA on mouse genes only
	genes <- rownames(object@raw.data)
} else { # mode == "hvg"
	# Compute the PCA on highly variable mouse genes
	if (length(object@var.genes) > 0) {
		genes <- object@var.genes
	} else {
		cat("Highly variable genes must be computed before PCA\n")
		q()
	}
}

object <- PCA(object = object, out.dir = OUT_DIR, pcs.compute = num_pc, suffix = mode, only_genes = genes, do.jack = jack)

save(object, file = paste(OUT_DIR, "/object.Robj", sep=""))


q()


