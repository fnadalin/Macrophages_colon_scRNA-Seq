# Global parameters

# resolution parameter for clustering
res <- 0.1

# number of PCs used
dim <- 0

# parameters for DEGs selection
min.pct <- 0.1
logFC <- 0.25

# parameters for DEGs filtering
logFC.filt <- 0.5
adjpval.filt <- 0.05
min.pct.cl <- 0.1

# number of top DEGs after filtering (for heatmap and volcano plots)
num_filt <- 50

filtered_cells <- NULL
cl_filter <- ""
renaming <- ""

width <- 5
height <- 7 

genes_highlight <- NULL

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
	cat("\nUsage: Rscript workflow_gene_expression.R <object> <out_dir> <mode> [<dim> <res> <cl_filter> <renaming> <width> <height> <genes_highlight>]\n\n")
	cat("<object>        .Robj file containing Seurat object \"obj_norm\"\n")
	cat("<out_dir>       where output will be stored\n")
	cat("<mode>          for clustering:\n")
	cat("                \"pc_mouse\": PCs computed on all mouse genes\n")
	cat("                \"hvg\": highly variable mouse genes\n")
	cat("                \"pc_mouse_hvg\": PCs computed on highly variable, non HIV-1 genes\n")
	cat("<dim>           number of PC to consider (only used in \"pc_notHIV\" or \"pc_hvg\" mode)\n")
	cat("<res>           resolution parameter for cluster identification\n")
	cat("<cl_filter>     comma-separated list of cluster IDs to be filtered out (OPTIONAL)\n")
	cat("<renaming>      comma-separated IDs to be assigned to the SELECTED clusters, after filtering (OPTIONAL)\n")
	cat("                Example: \"5,1,4,2,6,3\" means 0=>5, 1=>1, 2=>4, 3=>2, 4=>6, 5=>3\n")
	cat("<width>         width for the heatmap pdf plot\n")
	cat("<height>        height for the heatmap pdf plot\n")
	cat("<genes_highlight>  gene labels to be highlighted in the volcano pdf plot\n")
	q()
}

OBJECT <- args[1]
OUT_DIR <- args[2]
mode <- args[3]
if (length(args) > 3) {
	dim <- as.numeric(args[4])
	if (length(args) > 4) {
		res <- as.numeric(args[5])
		if (length(args) > 5) {
			cl_filter <- args[6]
			if (length(args) > 6) {
				renaming <- args[7]
				if (length(args) > 7) {
					width <- as.numeric(args[8])
					if (length(args) > 8) {
						height <- as.numeric(args[9])
						if (length(args) > 9) {
							GENES_HIGHLIGHT <- args[10]
							genes_highlight <- c(as.matrix(read.table(GENES_HIGHLIGHT)))
						}
					}
				}
			}
		}
	}
}
if (!(mode == "pc_mouse" || mode == "hvg" || mode == "pc_mouse_hvg")) {
	cat("<mode> must be either \"pc_mouse\" or \"hvg\" or \"pc_mouse_hvg\"\n")
	q()
}

if (mode == "pc_mouse" || dim == "pc_mouse_hvg") {
	if (length(args) < 4) {
		cat(paste("<dim> parameter is required for mode ", mode, "\n", sep=""))
		q()
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


object <- LoadObject(OBJECT)

prefix <- mode 
if (dim > 0) {
	prefix <- paste(prefix, dim, sep=".") 
}
id <- paste(prefix, "res", res, sep=".") 

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


# Find cluster gene markers
# tests to be used with normalized counts:
# - Wilcoxon (test whether two independent sets with the same variance have the same mean or not)
# - bimod (likelihood-ratio test for single cell gene expression)
# - MAST (GLM-framework that treates cellular detection rate as a covariate)
# tests to be used with raw counts:
# - poisson
# - negbinom
# - DESeq2

tests.norm <- c("wilcox", "MAST", "bimod")
tests.raw <- c("poisson", "negbinom") #, "DESeq2")
tests <- c(tests.norm, tests.raw)
tests <- c("MAST")


############# print the parameters to file ##############

cat("Print the parameters to file\n")

for (test.use in tests) {
	dir <- paste(OUT_DIR, "/", test.use, sep="")
	dir.create(dir, showWarnings = FALSE)
	sink(file = paste(dir, "/parameters.txt", sep = ""))
	cat(paste("res=", res, "\n", sep=""))
	cat(paste("dim=", dim, "\n", sep=""))
	cat(paste("min.pct=", min.pct, "\n", sep=""))
	cat(paste("logFC=", logFC, "\n", sep=""))
	cat(paste("logFC.filt=", logFC.filt, "\n", sep=""))
	cat(paste("adjpval.filt=", adjpval.filt, "\n", sep=""))
	cat(paste("min.pct.cl=", min.pct.cl, "\n", sep=""))
	cat(paste("num_filt=", num_filt, "\n", sep=""))
	cat(paste("cl_filter=\"", cl_filter, "\"\n", sep=""))
	cat(paste("renaming=\"", renaming, "\"\n", sep=""))
	sink()
}


######################## find DEGs ######################

cat("Find DEGs\n")

# cluster pairs
for (test.use in tests) {
	dir <- paste(OUT_DIR, "/", test.use, sep="")
	dir.create(dir, showWarnings = FALSE)
#	ClusterGeneMarkersByPairsNEW(object = object, out.dir = dir, res = res, prefix = prefix, clusters = new.cl.names, test.use = test.use, min.pct = min.pct, logFC = logFC)
#	ClusterGeneMarkersVsAllNEW(object, out.dir = dir, res = res, prefix = prefix,  clusters = new.cl.names, test.use = test.use, min.pct = min.pct, logFC = logFC)
}

####################### filter DEGs #####################

cat("Filter DEGs\n")

# filter out genes from the list to be plotted to the volcano plot
genes.use <- rownames(object@raw.data)[!(rownames(object@raw.data) %in% HIV_genes)]

# Filtering & heatmaps
# Output top 100 genes at most (top genes are defined according to the adj. p-value, hence they are method-specific)
# Volcano plot for all DEGs and label only the ones that pass the filtering step
for (test.use in tests) {

	# filter the DEGs based on FC and adj. p-value, and select the most significant 100 DEGs
	FilterClusterGeneMarkersPairsNEW(object, out.dir = OUT_DIR, res = res, prefix = prefix, clusters = new.cl.names, test.use = test.use, 
		logFC.filt = logFC.filt, adjpval.filt = adjpval.filt, min.pct1 = min.pct.cl, min.pct2 = min.pct.cl, num = num_filt)
#	HeatmapPlotPairsNEW(object, out.dir = OUT_DIR, test = test.use, clusters = new.cl.names, res = res, prefix = prefix)
	VolcanoPlotPairsNEW(object, out.dir = OUT_DIR, test = test.use, clusters = new.cl.names, res, prefix = prefix, genes.use = genes.use, highlight = genes_highlight)

	FilterClusterGeneMarkersAllNEW(object, out.dir = OUT_DIR, res = res, prefix = prefix, clusters = new.cl.names, test.use = test.use, 
		logFC.filt = logFC.filt, adjpval.filt = adjpval.filt, min.pct1 = min.pct.cl, min.pct2 = min.pct.cl, num = num_filt)
#	HeatmapPlotAllNEW(object, out.dir = OUT_DIR, test = test.use, clusters = new.cl.names, res = res, prefix = prefix, width = width, height = height)
	VolcanoPlotAllNEW(object, out.dir = OUT_DIR, test = test.use, clusters = new.cl.names, res, prefix = prefix, genes.use = genes.use, highlight = genes_highlight)

#	HeatmapPlotAllGlobal(object, out.dir = OUT_DIR, clusters = new.cl.names, test = test.use, res = res, prefix = prefix)
}


###################### DEG comparison ###################

q()

cat("DEG comparison (Pearson's correlation, Fisher exact test)\n")

for (t in c("tests.norm", "tests.raw")) {
	set <- eval(parse(text=t))

	# Compute the correlation between the adj. p-values across all DE methods
	DEGCorrelationPairs(out.dir = OUT_DIR, tests = set, clusters = new.cl.names, res = res, tests.name = t, plot = TRUE)
	DEGCorrelationAll(out.dir = OUT_DIR, tests = set, clusters = new.cl.names, res = res, tests.name = t, plot = TRUE)

	# Compute the significance of the intersection between the top significant DEGs across all DE methods
	DEGIntersectionPairs(out.dir = OUT_DIR, tests = set, clusters = new.cl.names, res = res, tests.name = t, plot = TRUE)
	DEGIntersectionAll(out.dir = OUT_DIR, tests = set, clusters = new.cl.names, res = res, tests.name = t, plot = TRUE)
}

q()


