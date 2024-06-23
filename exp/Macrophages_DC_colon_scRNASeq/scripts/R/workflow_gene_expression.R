# Global parameters

# mouse mito and ribo from:
# https://www.broadinstitute.org/files/shared/metabolism/mitocarta/mouse.mitocarta.2.0.html
# http://ribosome.med.miyazaki-u.ac.jp/

library("funr")

WORKING_DIR <- getwd()
SCRIPT_PATH <- dirname(sys.script())
SCRIPT_NAME <- basename(sys.script())
setwd(SCRIPT_PATH)
DATA_DIR <- "../../../../data"
setwd(WORKING_DIR)

MITO_GENES <- paste(DATA_DIR, "mouse_mitochondrial_genes.txt", sep="/")
RIBO_GENES <- paste(DATA_DIR, "mouse_ribosomal_genes.txt", sep="/")

dim <- 0
res <- 0.1

# seed for tSNE
seed <- 1

# perplexity for tSNE
perplexity <- 5

# parameters for UMAP
neighbors <- 30
min_dist <- 0.3

filtered_cells <- NULL
cl_filter <- ""
renaming <- ""

pt.size <- 1
label.size <- 12 # for clusters labelling on the tSNE/UMAP

# Signatures to be plotted cluster-wise

meta_genes <- c("nGene", "nUMI", "perc.mito.UMI", "perc.ribo.UMI")

# signatures <- c("CD11B", "CD11C", "CD103", "Dectin-1")
colon_signatures <- c("Itgam", "Itgax", "Itgae") # official gene symbols
epithelial <- c("Acta2", "Dcn", "Des", "Lamb1", "Lama4", "Lamc1", "Tpm2", "Nid1")
fibroblasts <- c("Car1", "Muc2", "Car4", "Krt19", "Krt8", "Krt18")
pathway <- c("Slc9a3r2", "Fxyd2", "Scnn1a", "Pik3r1", "Atp1b3", "Igf1", "Atp1a3")

violin_genes <- c("meta_genes", "colon_signatures", "epithelial", "fibroblasts", "pathway")
tsne_genes <- c("meta_genes", "colon_signatures", "epithelial", "fibroblasts", "pathway")


args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
	cat("\nUsage: Rscript workflow_gene_expression.R <object> <out_dir> <mode> [<dim> <res> <cl_filter> <renaming> <seed> <perplexity> <neighbors> <min_dist>]\n\n")
	cat("<object>        .Robj file containing Seurat object \"obj_norm\"\n")
	cat("<out_dir>       where output will be stored\n")
	cat("<mode>          for clustering:\n")
	cat("                \"pc_mouse\": PCs computed on all mouse genes\n")
	cat("                \"hvg\": highly variable mouse genes\n")
	cat("                \"pc_mouse_hvg\": PCs computed on highly variable, non HIV-1 genes\n")
	cat("<dim>           number of PC to consider (only used in \"pc_notHIV\" or \"pc_hvg\" mode)\n", sep="")
	cat("<res>           resolution parameter for cluster identification\n")
	cat("<cl_filter>     comma-separated list of cluster IDs to be filtered out (OPTIONAL)\n")
	cat("<renaming>      comma-separated IDs to be assigned to the SELECTED clusters, after filtering (OPTIONAL)\n")
	cat("                Example: \"5,1,4,2,6,3\" means 0=>5, 1=>1, 2=>4, 3=>2, 4=>6, 5=>3\n")
	cat(paste("<seed>          seed for the tSNE and the UMAP (default: ", seed, ")\n", sep=""))
	cat(paste("<perplexity>    tSNE perplexity (default: ", perplexity, ")\n", sep=""))
	cat(paste("<neighbors>     neighbors for the UMAP (default: ", neighbors, ")\n", sep=""))
	cat(paste("<min_dist>      min_dist for the UMAP (default: ", min_dist, ")\n\n", sep=""))
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
					seed <- as.numeric(args[8])
					if (length(args) > 8) {
						perplexity <- as.numeric(args[9])
						if (length(args) > 9) {
							neighbors <- as.numeric(args[10])
							if (length(args) > 10) {
								min_dist <- as.numeric(args[11])
							}
						}
					}
				}
			}
		}
	}
}

WORKING_DIR <- getwd()
SCRIPT_PATH <- dirname(sys.script())
SCRIPT_NAME <- basename(sys.script())
setwd(SCRIPT_PATH)
FUNCTIONS <- "functions.R"
source(FUNCTIONS)
setwd(WORKING_DIR)


# print the parameters to file
sink(file = paste(OUT_DIR, "/parameters_gene_expr.txt", sep = ""))
cat(paste("OBJECT=\"", OBJECT, "\"\n", sep=""))
cat(paste("OUT_DIR=\"", OUT_DIR, "\"\n", sep=""))
cat(paste("mode=\"", mode, "\"\n", sep=""))
cat(paste("dim=", dim, "\n", sep=""))
cat(paste("res=", res, "\n", sep=""))
cat(paste("cl_filter=\"", cl_filter, "\"\n", sep=""))
cat(paste("renaming=\"", renaming, "\"\n", sep=""))
cat(paste("seed=", seed, "\n", sep=""))
cat(paste("perplexity=", perplexity, "\n", sep=""))
cat(paste("neighbors=", neighbors, "\n", sep=""))
cat(paste("min_dist=", min_dist, "\n", sep=""))
cat(paste("pt.size=", pt.size, "\n", sep=""))
cat(paste("label.size=", label.size, "\n", sep=""))
sink()


object <- LoadObject(OBJECT)

prefix <- mode 
if (dim > 0) {
	prefix <- paste(prefix, dim, sep=".") 
}
tsne_reduction <- paste("tsne.", prefix, ".seed", seed, ".perp", perplexity, sep="")	
umap_reduction <- paste("umap.", prefix, ".seed", seed, ".neigh", neighbors, ".dist", min_dist, sep="")	

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

# generate the tSNE plots
cat("Generate the tSNE plots\n")
tSNEplotSamples (object = object, out.dir = OUT_DIR, title = "", seed = seed, pt.size = pt.size, reduction.use = tsne_reduction, cells.use = filtered_cells) 
tSNEplotClusters(object = object, out.dir = OUT_DIR, title = "", seed = seed, res = c(res), prefix = prefix, pt.size = pt.size, reduction.use = tsne_reduction, cells.use = filtered_cells, label.size = label.size)

# generate the UMAP plots
cat("Generate the UMAP plots\n")
tSNEplotSamples (object = object, out.dir = OUT_DIR, title = "", seed = seed, pt.size = pt.size, reduction.use = umap_reduction, cells.use = filtered_cells) 
tSNEplotClusters(object = object, out.dir = OUT_DIR, title = "", seed = seed, res = c(res), prefix = prefix, pt.size = pt.size, reduction.use = umap_reduction, cells.use = filtered_cells, label.size = label.size)

# clusters distribution on samples
SamplesComposition(object = object, out.dir = OUT_DIR, res = res, cl.names = new.cl.names, prefix = prefix)
# plot cluster distribution
BarplotSampleComposition(out.dir = OUT_DIR)

# GENE EXPRESSION

cat("Compute ribosomal and mitochondrial gene expression\n")
mito_genes <- c(as.matrix(read.table(MITO_GENES)))
ribo_genes <- c(as.matrix(read.table(RIBO_GENES)))
object <- GeneGroupExpression(object, genes = mito_genes, group_name = "mito", convert.to.official.name = TRUE, org="Mm")
object <- GeneGroupExpression(object, genes = ribo_genes, group_name = "ribo", convert.to.official.name = TRUE, org="Mm")
# object <- ComputeMeta(object = object)

cat("Plot gene expression (violin plots)\n")

# Plot gene expression across samples and clusters
for (name in violin_genes) {
	cat(paste(name, "\n"))
	genes <- eval(parse(text=name))
	PLOT_NAME <- paste(OUT_DIR, "/violin_samples_", name, ".pdf", sep="")
	ViolinPlotClusters(object = object, plot.name = PLOT_NAME, markers = genes, title = name, check = (length(grep(pattern="meta", x=name))==0))
	PLOT_NAME <- paste(OUT_DIR, "/violin_clusters_", name, ".pdf", sep="")
	ViolinPlotClusters(object = object, plot.name = PLOT_NAME, markers = genes, res = res, prefix = prefix, clusters = new.cl.names, title = name, check = (length(grep(pattern="meta", x=name))==0))
}

# Plot hvg
top20hvg <- head(rownames(object@hvg.info),20)
PLOT_NAME <- paste(OUT_DIR, "/violin_samples_top20hvg.pdf", sep="")
ViolinPlotClusters(object, plot.name = PLOT_NAME, markers = top20hvg, res = res, prefix = prefix, clusters = new.cl.names, title = "top 20 hvg")

cat("Plot gene expression (tSNE/UMAP plots)\n")

colors <- c("gray", "red")

for (name in tsne_genes) {
	cat(paste(name, "\n"))
	genes <- eval(parse(text=name))
	for (gene in genes) {
		if (gene %in% c(object@data@Dimnames[[1]], colnames(object@meta.data))) {
			for (reduction in c(tsne_reduction, umap_reduction)) {
				PLOT_NAME <- paste(OUT_DIR, "/", reduction, "_", gene, ".pdf", sep="")
				pdf(PLOT_NAME, width=6, height=5)
				FeaturePlot(object = object, features.plot = gene, reduction.use = reduction, pt.size = pt.size, no.axes=TRUE, no.legend = FALSE, cells.use = filtered_cells, cols.use = colors)
				dev.off()
			}
		}
	}
}

q()


