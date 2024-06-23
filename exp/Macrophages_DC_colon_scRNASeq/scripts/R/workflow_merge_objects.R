# Select the row data from a cluster on the first object
# Add the second sample (apply the constraints for ngenes and ncells)
# Normalize and scale

library("funr")

WORKING_DIR <- getwd()
SCRIPT_PATH <- dirname(sys.script())
SCRIPT_NAME <- basename(sys.script())
setwd(SCRIPT_PATH)
DATA_DIR <- "../..data"
setwd(WORKING_DIR)

MITO_GENES <- paste(DATA_DIR, "mouse_mitochondrial_genes.txt", sep="/")
RIBO_GENES <- paste(DATA_DIR, "mouse_ribosomal_genes.txt", sep="/")

mito_genes <- c(as.matrix(read.table(MITO_GENES)))
ribo_genes <- c(as.matrix(read.table(RIBO_GENES)))

ORG <- "Mm"

min_cells <- 3
min_genes <- 200

# limits for the average and dispersion for the detection of highly variable genes
xmin <- 0.1
xmax <- 8
ymin <- 1
# xmin = 0.0125, xmax = 3, ymin = 0.5 (Seurat tutorial)

violin_plot_features <- c("nUMI", "nGene", "perc.mito.UMI", "perc.ribo.UMI")

############################ IMPORTANT #############################
#   Seurat requires cells >= min_cells and genes > min_genes!!!!   #
####################################################################

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
	cat("\nUsage: Rscript workflow_merge_objects.R <object1> <object2> <out_dir> [<min_cells> <min_genes> <xmin> <xmax> <ymin>]\n")
	cat("<object1>        object 1 file name\n")
	cat("<object2>        object 2 file name\n")
	cat(paste("<ncell>         minimum number of cells per gene [OPTIONAL - default: ", min_cells, "]\n", sep=""))
	cat(paste("<ngene>         minimum number of genes per cell [OPTIONAL - default: ", min_genes, "]\n", sep=""))
	cat(paste("<xmin>          min average expression for highly variable genes [OPTIONAL - default: ", xmin, "]\n", sep=""))
	cat(paste("<xmax>          max average expression for highly variable genes [OPTIONAL - default: ", xmax, "]\n", sep=""))
	cat(paste("<ymin>          min dispersion for highly variable genes [OPTIONAL - default: ", ymin, "]\n\n", sep=""))
	q()
}

OBJECT1 <- args[1]
OBJECT2 <- args[2]
OUT_DIR <- args[3]
if (length(args) > 3) {
	min_cells <- as.numeric(args[4])
	if (length(args) > 4) {
		min_genes <- as.numeric(args[5])
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
sink(file = paste(OUT_DIR, "/parameters_merge_objects.txt", sep = ""))
cat(paste("OBJECT1=", OBJECT1, "\n", sep=""))
cat(paste("OBJECT2=", OBJECT2, "\n", sep=""))
cat(paste("OUT_DIR=", OUT_DIR, "\n", sep=""))
cat(paste("min_cells=", min_cells, "\n", sep=""))
cat(paste("min_genes=", min_genes, "\n", sep=""))
cat(paste("max_genes=", min_genes, "\n", sep=""))
cat(paste("perc.mito=", perc.mito, "\n", sep=""))
cat(paste("xmin=", xmin, "\n", sep=""))
cat(paste("xmax=", xmax, "\n", sep=""))
cat(paste("ymin=", ymin, "\n", sep=""))
sink()

object1 <- LoadObject(OBJECT1)
object2 <- LoadObject(OBJECT2)

# add the sample
ident1 <- object1@ident
ident2 <- object2@ident
object <- MergeSeurat(object1 = object1, object2 = object2, add.cell.id1 = ident1, add.cell.id2 = ident2, project="DC", do.normalize=FALSE, min.cells = min_cells, min.genes = (min_genes-1))

# Compute number of mito and ribo genes
object <- GeneGroupExpression(object, genes = mito_genes, group_name = "mito", convert.to.official.name = TRUE, org=ORG)
object <- GeneGroupExpression(object, genes = ribo_genes, group_name = "ribo", convert.to.official.name = TRUE, org=ORG)

# Plot meta data after filtering
PlotMetaData(object, violin_plot_features, OUT_DIR, filtered = TRUE)

# normalize and scale
object <- NormScaleMatrix(object = object, model="linear", vars.to.regress="nUMI")

object <- FindVariableGenes(object, x.low.cutoff = xmin, x.high.cutoff = xmax, y.cutoff = ymin)

save(object, file = paste(OUT_DIR, "/object.Robj", sep=""))




q()

