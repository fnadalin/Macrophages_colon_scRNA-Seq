# Global parameters

res <- NULL


args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
	cat("\nUsage: Rscript workflow_sample_composition.R <object> <out_dir> <mode> <res> [<dim>]\n\n")
	cat("<object>        .Robj file containing Seurat object \"obj_norm\"\n")
	cat("<out_dir>       where output will be stored\n")
	cat("<mode>          for clustering:\n")
	cat("                \"pc_mouse\": PCs computed on all mouse genes\n")
	cat("                \"hvg\": highly variable mouse genes\n")
	cat("                \"pc_mouse_hvg\": PCs computed on highly variable mouse genes\n")
	cat("<res>           resolution parameter for clustering\n")
	cat(paste("<dim>           number of PC to consider (only used in \"pc_notHIV\" or \"pc_hvg\" mode)\n\n", sep=""))
	q()
}

OBJECT <- args[1]
OUT_DIR <- args[2]
mode <- args[3]
res <- eval(parse(text=args[4]))

if (!(mode == "pc_mouse" || mode == "hvg" || mode == "pc_mouse_hvg")) {
	cat("<mode> must be either \"pc_mouse\" or \"hvg\" or \"pc_mouse_hvg\"\n")
	q()
}

if (mode == "pc_mouse" || mode == "pc_mouse_hvg") {
	if (length(args) < 5) {
		cat(paste("<dim> parameter is required for mode ", mode, "\n", sep=""))
		q()
	} else {
		dim <- as.numeric(args[5])
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
sink(file = paste(OUT_DIR, "/parameters_samples_composition.txt", sep = ""))
cat(paste("object=", OBJECT, "\n", sep=""))
cat(paste("out.dir=", OUT_DIR, "\n", sep=""))
cat(paste("mode=", mode, "\n", sep=""))
cat(paste("res=", res, "\n", sep=""))
cat(paste("dim=", dim, "\n", sep=""))
sink()


object <- LoadObject(OBJECT)

prefix <- mode 
if (dim > 0) {
	prefix <- paste(prefix, dim, sep=".") 
}
reduction <- paste("samples.composition.", prefix, sep="")

out_dir_samples <- paste(OUT_DIR, reduction, sep="/")
dir.create(out_dir_samples, showWarnings = FALSE)

for (r in res) {
	SamplesComposition(object, out.dir = out_dir_samples, res = r, prefix = prefix)
}

BarplotSampleCompositionMulti(object, out.dir = out_dir_samples, res, prefix = prefix) 


q()


