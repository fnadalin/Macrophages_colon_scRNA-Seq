# Created on 2019/10/10

########################### DEFAULT PARAMETER VALUES ###########################

# for cell labelling
ADD_ID <- FALSE

# for gene/cell filtering
MIN_CELLS <- 3
MIN_CELLS_PER_SAMPLE <- 0
MIN_GENES <- 200
MAX_GENES <- 100000
PERC_MITO <- 0.0

# limits for the average and dispersion for the detection of highly variable genes
XMIN <- 0.1
XMAX <- 8
YMIN <- 1
# xmin = 0.0125, xmax = 3, ymin = 0.5 (Seurat tutorial)

RegressUMI <- FALSE

DOUBLET_LIST <- ""

GENE_GROUP_FILE <- NULL
GENE_GROUP_ID <- NULL

outliers <- FALSE
STD_FACTOR <- 2
bimod <- FALSE

CONTAMINANT_FILE <- NULL
CONTAMINANT_ID <- NULL

CONTAMINANT_STD_FACTOR <- 2
use_auc <- FALSE
AUC_MAX_RANK_PERC <- 0.2
AUC_CUTOFF <- 0.5

violin_plot_meta <- c("nCount_RNA", "nFeature_RNA")
VIOLIN_PT_SIZE <- 1


################################# OPTION MENU ##################################

library("optparse")

# specify our desired options in a list
# by default OptionParser will add an help option equivalent to 
# make_option(c("-h", "--help"), action="store_true", default=FALSE, 
#               help="Show this help message and exit")
cat("\n")
option_list <- list( 
    make_option("--object1", type = "character",
        help="[REQUIRED] .Robj file containing the first Seurat object"),
    make_option("--object2", type = "character",
        help="[REQUIRED] .Robj file containing the second Seurat object"),
    make_option("--out_dir", type = "character",
        help="[REQUIRED] directory where output will be stored"),
	make_option("--add_id", action = "store_true",
		help="paste the identity of the single objects at the beginning of barcode IDs (only if the active identity of the cells in the same object is unique!)"),
    make_option("--min_cells", type = "integer", default = MIN_CELLS,
        help="minimum number of cells per gene [default=%default]"),
    make_option("--min_cells_per_sample", type = "integer", default = MIN_CELLS_PER_SAMPLE,
        help="minimum number of cells per gene, for each sample [retrieved from object@meta.data$orig.ident - default=%default]"),
    make_option("--min_genes", type = "integer", default = MIN_GENES,
        help="minimum number of genes per cell [default=%default]"),
    make_option("--max_genes", type = "integer", default = MAX_GENES,
        help="maximum number of genes per cell [default=%default]"),
	make_option("--outliers", action = "store_true",
        help="remove cells where the number of detected genes is more than STD_FACTOR standard deviations away from the mean"),
	make_option("--std_factor", type = "double", default = STD_FACTOR,
        help="number of standard deviations to identify outliers [default=%default]"),
	make_option("--bimod", action = "store_true",
        help="remove cells with lower gene count in a bimodal distribution (can be used in combination with --outliers)"),
	make_option("--doublet_list", type = "character",
        help="names of the files with values indicating predicted doublets for each sample, one per line"),
    make_option("--gene_group_files", type = "character",
        help="comma-separated list of file names containing the gene symbols, one per line, which are used to compute a joint expression level"),
    make_option("--gene_group_IDs", type = "character",
        help="comma-separated list of gene group IDs"),
	make_option("--perc_mito", type = "double", default = PERC_MITO, 
        help="maximum fraction of mitochondrial UMIs per cell (requires the computation of mitochondrial gene group expression \"mito.UMI\") [default=%default]"),
    make_option("--contaminant_files", type = "character",
        help="comma-separated list of file names containing the gene markers, one per line, for each contaminant cell type"),
    make_option("--contaminant_IDs", type = "character",
        help="comma-separated list of contaminant cell type IDs"),
	make_option("--contaminant_std_factor", type = "double", default = CONTAMINANT_STD_FACTOR,
        help="remove cells where the expression of any contaminant gene is more than STD_FACTOR standard deviations above the mean [default=%default]"),
	make_option("--use_auc", action = "store_true",
        help="use AUCell to filter out contaminants (overrides contaminant filtering based on outliers identification)"),	
	make_option("--auc_max_rank_perc", type = "double", default = AUC_MAX_RANK_PERC, 
        help="used to compute aucMaxRank value in AUCell:::AUCell_calcAUC() [default=%default]"),
	make_option("--auc_cutoff", type = "double", default = AUC_CUTOFF, 
        help="high.thresholds value in Seurat:::FilterCells() [default=%default]"),	
    make_option("--violin_pt_size", type = "double", default = VIOLIN_PT_SIZE,
        help="point size for violin plots [default=%default]"),
    make_option("--regressUMI", action = "store_true",
        help="regress out the UMI count before scaling"),
    make_option("--xmin", type = "double", default = XMIN,
        help="min average expression for highly variable genes [default=%default]"),
    make_option("--xmax", type = "double", default = XMAX,
        help="max average expression for highly variable genes [default=%default]"),
    make_option("--ymin", type = "double", default = YMIN,
        help="min dispersion for highly variable genes [default=%default]")
)


################################ PARSE OPTIONS #################################

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults
opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$object1)) {
	write("Option --object1 is required\nTry --help for help", stderr()) 
	q()
} else {
	OBJECT1 <- opt$object1
}

if (is.null(opt$object2)) {
	write("Option --object2 is required\nTry --help for help", stderr()) 
	q()
} else {
	OBJECT2 <- opt$object2
}

if (is.null(opt$out_dir)) {
	write("Option --out_dir is required\nTry --help for help", stderr()) 
	q()
} else {
	OUT_DIR <- opt$out_dir
}

if (!is.null(opt$add_id)) {
	ADD_ID <- opt$add_id
}

if (!is.null(opt$min_cells))
	MIN_CELLS <- opt$min_cells

if (!is.null(opt$min_cells_per_sample))
	MIN_CELLS_PER_SAMPLE <- opt$min_cells_per_sample

if (!is.null(opt$min_genes))
	MIN_GENES <- opt$min_genes

if (!is.null(opt$max_genes))
	MAX_GENES <- opt$max_genes

if (!is.null(opt$outliers))
	outliers <- TRUE

if (!is.null(opt$std_factor))
	STD_FACTOR <- opt$std_factor

if (!is.null(opt$bimod))
	bimod <- TRUE

if (!is.null(opt$doublet_list))
	DOUBLET_LIST <- opt$doublet_list

if (!is.null(opt$gene_group_files))
	GENE_GROUP_FILE <- unlist(strsplit(x = opt$gene_group_files, split = ","))

if (!is.null(opt$gene_group_IDs))
	GENE_GROUP_ID <- unlist(strsplit(x = opt$gene_group_IDs, split = ","))

if (!is.null(opt$perc_mito))
	PERC_MITO <- opt$perc_mito

if (!is.null(opt$contaminant_files))
	CONTAMINANT_FILE <- unlist(strsplit(x = opt$contaminant_files, split = ","))

if (!is.null(opt$contaminant_IDs))
	CONTAMINANT_ID <- unlist(strsplit(x = opt$contaminant_IDs, split = ","))

if (!is.null(opt$use_auc)) 
	use_auc <- TRUE

if (!is.null(opt$auc_max_rank_perc))
	AUC_MAX_RANK_PERC <- opt$auc_max_rank_perc

if (!is.null(opt$auc_cutoff))
	AUC_CUTOFF <- opt$auc_cutoff

if (!is.null(opt$regressUMI)) 
	RegressUMI <- TRUE

if (!is.null(opt$xmin)) 
	XMIN <- opt$xmin

if (!is.null(opt$xmax)) 
	XMAX <- opt$xmax

if (!is.null(opt$ymin)) 
	YMIN <- opt$ymin

if (!is.null(opt$violin_pt_size)) 
	VIOLIN_PT_SIZE <- opt$violin_pt_size

if (length(GENE_GROUP_FILE) != length(GENE_GROUP_ID)) {
	write("Options --geneGroupFiles and --geneGroupIDs must contain the same number of groups\nTry --help for help", stderr()) 
	q()
}


############################### EXPORT FUNCTIONS ###############################

library("funr")

WORKING_DIR <- getwd()
SCRIPT_PATH <- dirname(sys.script())
SCRIPT_NAME <- basename(sys.script())
setwd(SCRIPT_PATH)
FUNCTIONS <- "functions.R"
source(FUNCTIONS)
setwd(WORKING_DIR)

######################### PRINT THE PARAMETERS TO FILE #########################

dir.create(OUT_DIR, showWarnings = FALSE)
sink(file = paste(OUT_DIR, "/parameters.txt", sep = ""))
cat("================ Load parameters ================\n")
for (i in 1:length(opt))
	cat(paste(names(opt)[i], "=", opt[[i]], "\n", sep=""))
sink()


################################## EXECUTION ###################################


object1 <- LoadObject(OBJECT1)
object2 <- LoadObject(OBJECT2)

ident1 <- as.character(unique(Idents(object1)))
ident2 <- as.character(unique(Idents(object2)))

if (!ADD_ID || length(ident1) > 1 || length(ident2) > 1) {
	object <- merge(x = object1,
                    y = object2,
                    min.cells = MIN_CELLS,
                    min.features = (MIN_GENES-1))
} else {
	object <- merge(x = object1,
                    y = object2,
                    min.cells = MIN_CELLS,
                    min.features = (MIN_GENES-1),
                    add.cell.ids = c(ident1,ident2))
}

if (MIN_CELLS_PER_SAMPLE > 0) {
	id <- object$orig.ident
	id_list <- unique(id)
	v <- rep(0,	length(rownames(object)))
	counts <- GetAssayData(object, slot="counts")
	for (x in id_list) {
		num.cells <- rowSums(counts[,id==x] > 0)
		v <- v + (num.cells >= MIN_CELLS_PER_SAMPLE)
	}
	genes.use <- rownames(object)[v==length(id_list)]
	subset.matrix <- object@raw.data[genes.use, ]
	object2 <- CreateSeuratObject(subset.matrix) # Create a new Seurat object with just the genes of interest
	meta <- object@meta.data
	apply(1:ncol(meta), 2, function(x) { 
			v <- meta[,x]
			names(v) <- colnames(x = object) 
			object2 <- AddMetaData(x = object2, metadata = v, col.name = colnames(meta)[x]) }
	)
}

PlotMetaData(object, violin_plot_meta, OUT_DIR, point_size = VIOLIN_PT_SIZE)

if (!is.null(GENE_GROUP_FILE)) {
	for (i in 1:length(GENE_GROUP_FILE)) {
		genes <- as.matrix(read.table(GENE_GROUP_FILE[i]))
		object <- GeneGroupExpression(object = object, 
		                              genes = genes, 
		                              group_name = GENE_GROUP_ID[i])
	}
	violin_plot_features <- gsub("$", ".UMI", gsub("^", "perc.", GENE_GROUP_ID))
	PlotMetaData(object, violin_plot_features, OUT_DIR, point_size = VIOLIN_PT_SIZE)
}

object <- BimodalnGeneFilter(object = object, 
                             bimod = bimod, 
                             outliers = outliers, 
                             perc.mito = PERC_MITO, 
                             min.genes = MIN_GENES,
                             max.genes = MAX_GENES, 
                             second_max_interval_start = SECOND_MAX_INTERVAL_START,
                             second_max_interval_end_subtract = SECOND_MAX_INTERVAL_END_SUBTRACT,
                             sd_factor = STD_FACTOR)

# Normalize UMI counts by cell 
# Regress out the number of UMIs per cell
# Center and scale
if (RegressUMI) {
	object <- NormScaleMatrix(object = object)
} else {
	object <- NormScaleMatrix(object = object, 
                              vars.to.regress = NULL)
}

# now the normalized expression has been computed
if (!is.null(CONTAMINANT_FILE)) {
	cat("Filter by contaminants\n")
	if (!use_auc) {
		object <- FilterCellsByHighGeneExpression(object = object, 
                                                  genes = CONTAMINANT_FILE, 
                                                  sd_factor = CONTAMINANT_STD_FACTOR)
	} else {
		geneSets <- as.list(apply(CONTAMINANT_FILE, function(x) c(as.matrix(read.table(x)))))
		names(geneSets) <- CONTAMINANT_ID
		object <- FilterCellsByGeneAuc(object = object, 
                                       genes = geneSets, 
                                       aucMaxRankPerc = AUC_MAX_RANK_PERC, 
                                       auc.cutoff = AUC_CUTOFF)
	}
}

PlotMetaData(object, violin_plot_meta, OUT_DIR, filtered = TRUE, point_size = VIOLIN_PT_SIZE)

if (!is.null(GENE_GROUP_FILE)) {
	violin_plot_features <- gsub("$", ".UMI", gsub("^", "perc.", GENE_GROUP_ID))
	PlotMetaData(object, violin_plot_features, OUT_DIR, filtered = TRUE, point_size = VIOLIN_PT_SIZE)
}

object <- FindVariableFeatures(object = object, 
                               mean.cutoff = c(XMIN, XMAX), 
                               dispersion.cutoff = c(YMIN, Inf))

save(object, file = paste(OUT_DIR, "/object.Robj", sep=""))

sessionInfo()


q()

