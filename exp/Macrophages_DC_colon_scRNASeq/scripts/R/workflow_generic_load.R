
library("funr")

WORKING_DIR <- getwd()
SCRIPT_PATH <- dirname(sys.script())
SCRIPT_NAME <- basename(sys.script())
setwd(SCRIPT_PATH)
DATA_DIR <- "../../../../data"
setwd(WORKING_DIR)

MITO_GENES <- paste(DATA_DIR, "mouse_mitochondrial_genes.txt", sep="/")
RIBO_GENES <- paste(DATA_DIR, "mouse_ribosomal_genes.txt", sep="/")

mito_genes <- c(as.matrix(read.table(MITO_GENES)))
ribo_genes <- c(as.matrix(read.table(RIBO_GENES)))

ORG <- "Mm" # for gene name conversion 

min_cells <- 3
min_genes <- 200
max_genes <- 10000
perc.mito <- 0.1

# limits for the average and dispersion for the detection of highly variable genes
xmin <- 0.1
xmax <- 8
ymin <- 1
# xmin = 0.0125, xmax = 3, ymin = 0.5 (Seurat tutorial)

outliers <- FALSE
bimod <- FALSE
second_max_interval_start <- 200
second_max_interval_end_subtract <- 1

contaminants <- c()

STD_FACTOR <- 2
auc <- FALSE
AUC_MAX_RANK_PERC <- 0.2
AUC_CUTOFF <- 0.5

violin_plot_features <- c("nUMI", "nGene", "perc.mito.UMI", "perc.ribo.UMI")

# signatures <- c("CD11B", "CD11C", "CD103", "Dectin-1")
colon_signatures <- c("Itgam", "Itgax", "Itgae") # official gene symbols
epithelial <- c("Acta2", "Dcn", "Des", "Lamb1", "Lama4", "Lamc1", "Tpm2", "Nid1")
fibroblasts <- c("Car1", "Muc2", "Car4", "Krt19", "Krt8", "Krt18")
pathway <- c("Slc9a3r2", "Fxyd2", "Scnn1a", "Pik3r1", "Atp1b3", "Igf1", "Atp1a3")

contaminant_cells <- c("epithelial", "fibroblasts")
violin_genes <- c("colon_signatures", "epithelial", "fibroblasts", "pathway")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
	cat("Usage: Rscript workflow_generic_load.R <files_list> <names_list> <out_dir> [<ncell> <ngene> <xmin> <xmax> <ymin> <perc.mito> <outliers_nGene> <bimod_nGene> <second_max_interval_start> <second_max_interval_end_subtract> <contaminants>]\n")
	cat("<files_list>    name of the directories where CellRanger-count output is stored for each sample, one per line (e.g. out_Cellranger/H12-CpGA_multi/GRCh38/)\n")
	cat("<names_list>    names to assign to each sample (N.B: it is important that they do not contain \"_\")\n")
	cat("<out_dir>       where output will be stored\n")
	cat(paste("<ncell>            minimum number of cells per gene [OPTIONAL - default: ", min_cells, "]\n", sep=""))
	cat(paste("<ngene>            minimum number of genes per cell [OPTIONAL - default: ", min_genes, "]\n", sep=""))
	cat(paste("<xmin>             min average expression for highly variable genes [OPTIONAL - default: ", xmin, "]\n", sep=""))
	cat(paste("<xmax>             max average expression for highly variable genes [OPTIONAL - default: ", xmax, "]\n", sep=""))
	cat(paste("<ymin>             min dispersion for highly variable genes [OPTIONAL - default: ", ymin, "]\n", sep=""))
	cat(paste("<perc.mito>        max percentage of mitochondrial genes [OPTIONAL - default: ", perc.mito, "]\n", sep=""))
	cat(paste("<outliers_nGene>   remove cells with lower/higher nGene count compared to the distribution [OPTIONAL - default: no filtering]\n", sep=""))
	cat(paste("<bimod_nGene>      remove cells with lower nGene count in a bimodal (\"b\") distribution [OPTIONAL - default: no filtering]\n", sep=""))
	cat("<contaminants>     remove cells with higher expression (compared to the mean) for either one of the genes specified in the list\n")
	cat("<use_auc>          instead of discarding cells based on gene outliers, use AUCell\n\n")
	q()
}

FILES_LIST <- args[1]
NAMES_LIST <- args[2]
OUT_DIR <- args[3]
if (length(args) > 3) {
	min_cells <- as.numeric(args[4])
	if (length(args) > 4) {
		min_genes <- as.numeric(args[5])
		if (length(args) > 5) {
			xmin <- as.numeric(args[6])
			if (length(args) > 6) {
				xmax <- as.numeric(args[7])
				if (length(args) > 7) {
					ymin <- as.numeric(args[8])
					if (length(args) > 8) {
						perc.mito <- as.numeric(args[9])
						if (length(args) > 9) {
							outliers <- (args[10] == 'o')
							if (length(args) > 10) {
								bimod <- (args[11] == "b")
								if (length(args) > 11) {
									second_max_interval_start <- as.numeric(args[12])
									if (length(args) > 12) {
										second_max_interval_end_subtract <- as.numeric(args[13])
										if (length(args) > 13) {
											CONTAMINANTS <- args[14]
											contaminants <- c(as.matrix(read.table(CONTAMINANTS)))
											if (length(args) > 14) {
												auc <- (args[15] == "a")
											}
										}
									}
								}
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
sink(file = paste(OUT_DIR, "/parameters.txt", sep = ""))
cat(paste("FILES_LIST=", FILES_LIST, "\n", sep=""))
cat(paste("NAMES_LIST=", NAMES_LIST, "\n", sep=""))
cat(paste("OUT_DIR=", OUT_DIR, "\n", sep=""))
cat(paste("min_cells=", min_cells, "\n", sep=""))
cat(paste("min_genes=", min_genes, "\n", sep=""))
cat(paste("max_genes=", min_genes, "\n", sep=""))
cat(paste("perc.mito=", perc.mito, "\n", sep=""))
cat(paste("xmin=", xmin, "\n", sep=""))
cat(paste("xmax=", xmax, "\n", sep=""))
cat(paste("ymin=", ymin, "\n", sep=""))
cat(paste("outliers=", outliers, "\n", sep=""))
cat(paste("bimod=", bimod, "\n", sep=""))
cat(paste("second_max_interval_start=", second_max_interval_start, "\n", sep=""))
cat(paste("second_max_interval_end_subtract=", second_max_interval_end_subtract, "\n", sep=""))
cat(paste("STD_FACTOR=", STD_FACTOR, "\n", sep=""))
cat(paste("auc=", auc, "\n", sep=""))
cat(paste("AUC_MAX_RANK_PERC=", AUC_MAX_RANK_PERC, "\n", sep=""))
cat(paste("AUC_CUTOFF=", AUC_CUTOFF, "\n", sep=""))
sink()

# Load the UMI counts 
object <- LoadMultiMatrix(files.list = FILES_LIST, names.list = NAMES_LIST, out.dir = OUT_DIR, min.cells=min_cells, min.genes=(min_genes-1))

# Compute number of mito and ribo genes
object <- GeneGroupExpression(object, genes = mito_genes, group_name = "mito", convert.to.official.name = TRUE, org=ORG)
object <- GeneGroupExpression(object, genes = ribo_genes, group_name = "ribo", convert.to.official.name = TRUE, org=ORG)

# Plot meta data before filtering
PlotMetaData(object, violin_plot_features, OUT_DIR)

# filter by low/high gene count and high perc of mito
object <- BimodalnGeneFilter(object, bimod = bimod, outliers = TRUE, perc.mito = perc.mito, min.genes=min_genes, max.genes = max_genes, second_max_interval_start = second_max_interval_start, 
			second_max_interval_end_subtract = second_max_interval_end_subtract, sd_factor = STD_FACTOR)

# Plot meta data after filtering
PlotMetaData(object, violin_plot_features, OUT_DIR, filtered = TRUE)

# Normalize UMI counts by cell 
# Regress out the number of UMIs per cell
# Center and scale
object <- NormScaleMatrix(object)

# now the normalized expression has been computed
if (!is.null(contaminants)) {
	cat("Filter by contaminants\n")
	if (!auc) {
		object <- FilterCellsByHighGeneExpression(object, genes = contaminants, sd_factor = STD_FACTOR)
	} else {
		geneSets <- list()
		for (i in 1:length(contaminant_cells)) {
			geneSets[[i]] <- eval(parse(text=contaminant_cells[i]))
		}
		names(geneSets) <- contaminant_cells
		object <- FilterCellsByGeneAuc(object, genes = geneSets, aucMaxRankPerc = AUC_MAX_RANK_PERC, auc.cutoff = AUC_CUTOFF)
	}
}

PLOT_NAME <- paste(OUT_DIR, "hvg.pdf", sep="/")
pdf(PLOT_NAME)
object <- FindVariableGenes(object, x.low.cutoff = xmin, x.high.cutoff = xmax, y.cutoff = ymin)
dev.off()

save(object, file = paste(OUT_DIR, "/object.Robj", sep=""))


q()


