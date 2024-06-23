# Global parameters

# resolution parameter for clustering
res <- 0.1

# number of PCs used
dim <- 0

filtered_cells <- NULL
cl_filter <- ""
renaming <- ""

# criterion used for clustering
# might not be specified at all: in such a case, the corresponding field in Seurat
# object is assigned the default name, with no prefix
mode <- ""

prefix <- NULL

# measure to be used to rank the genes
measure <- "t" # "f" for fold change


args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
	cat("\nUsage: Rscript workflow_annotation_preprocessing.R <object> <out_dir> [<mode> <measure> <dim> <res> <cl_filter> <renaming>]\n\n")
	cat("<object>        .Robj file containing Seurat object \"obj_norm\"\n")
	cat("<out_dir>       where output will be stored\n")
	cat("<mode>          for clustering:\n")
	cat("                \"pc_mouse\": PCs computed on all mouse genes\n")
	cat("                \"hvg\": highly variable mouse genes\n")
	cat("                \"pc_mouse_hvg\": PCs computed on highly variable mouse genes\n")
	cat(paste("<measure>       measure to be used to rank genes [OPTIONAL - default: ", measure, "]\n", sep=""))
	cat("<dim>           number of PC to consider [OPTIONAL - only used in \"pc_notHIV\" or \"pc_hvg\" mode]\n")
	cat(paste("<res>           resolution parameter for cluster identification [OPTIONAL - default: ", res, "]\n", sep=""))
	cat("<cl_filter>     comma-separated list of cluster IDs to be filtered out [OPTIONAL]\n")
	cat("<renaming>      comma-separated IDs to be assigned to the SELECTED clusters, after filtering [OPTIONAL]\n")
	cat("                Example: \"5,1,4,2,6,3\" means 0=>5, 1=>1, 2=>4, 3=>2, 4=>6, 5=>3\n\n")
	q()
}

OBJECT <- args[1]
OUT_DIR <- args[2]
if (length(args) > 2) {
	mode <- args[3]
		if (length(args) > 3) {
		measure <- args[4]
		if (length(args) > 4) {
			dim <- as.numeric(args[5])
			if (length(args) > 5) {
				res <- as.numeric(args[6])
				if (length(args) > 6) {
					cl_filter <- args[7]
					if (length(args) > 7) {
						renaming <- args[8]
					}
				}
			}
		}
	}
}

if (mode != "") {
	if (!(mode == "pc_mouse" || mode == "hvg" || mode == "pc_mouse_hvg")) {
		cat("<mode> must be either \"pc_mouse\" or \"hvg\" or \"pc_mouse_hvg\"\n")
		q()
	}
	if (mode == "pc_mouse" || mode == "pc_mouse_hvg") {
		if (dim == 0) {
			cat(paste("<dim> parameter is required for mode ", mode, "\n", sep=""))
			q()
		}
	}
	prefix <- mode
	if (dim > 0) {
		prefix <- paste(prefix, dim, sep=".") 
	}
}

PATH="~/curie/exp/Macrophages_DC_colon_scRNASeq/scripts/R"
source(paste(PATH, "/functions.R", sep=""))


object <- LoadObject(OBJECT)

id <- ClustersID(res, prefix)

if (cl_filter != "") {
	# filter out cells
	filtered_cells <- FilterOutClusters (object, id, cl_filter)
}

# rename the clusters
object <- RenameClusters(object, id, renaming, cl_filter)
if (renaming == "") {
	new.cl.names <- ClusterNames(object, id, exclude = c("X"))
} else {
	new.cl.names <- unlist(strsplit(renaming, split=","))
}


###################### Sort genes according to != measures #####################

cat("Preprocessing for GSEA\n")

dir <- paste(OUT_DIR, "GSEAprep", sep="/")
dir.create(dir, showWarnings = FALSE)
GSEAPreprocessingClusterByPairs(object = object, out.dir = dir, res = res, prefix = prefix, clusters = new.cl.names, measure = measure)
GSEAPreprocessingClusterVsAll(object, out.dir = dir, res = res, prefix = prefix,  clusters = new.cl.names, measure = measure)


q()


