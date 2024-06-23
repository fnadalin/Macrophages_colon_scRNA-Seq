########## NB ###########
# In order to plot a pdf in a loop or from
# withi a function, the plot command must 
# be surrounded by the print() function!!!
#########################


library("Seurat")
library("RColorBrewer")
library("ggplot2")
library("grid")
library("ggrepel")
library("scales") # to emulate ggplot default color palette
library("factoextra") # for clusters validation
library("cluster") # for clusters validation
library("pbapply")
library("scater") # for MAST use in Seurat latest version
library("limma") # for alias2SymbolTable
library("AUCell") # for cells filtering


MINIMUM <- 1e-300
MAXIMUM <- 100000
EPSILON <- 0.0001

all_names <- c("H12-mock", "H24-mock", "H12-HIV", "H24-HIV", "H12-CpGA", "H24-CpGA", "H24-HIVAztNvp")

# How to order the sample names
# ordering <- c("H12-mock", "H12-HIV", "H24-mock", "H24-HIVAztNvp", "H24-HIV")
ordering <- c("H12-mock", "H12-HIV", "H12-CpGA", "H24-mock", "H24-HIVAztNvp", "H24-HIV", "H24-CpGA")

colors <- brewer.pal(8, "Paired")[c(1:6,8)]


######### groups of genes ##########

HIV_genes <- c("gag-pol", "vif", "vpr", "tat", "rev", "vpu", "env", "asp", "nef")

preDC_genes_for_gating <- preDC_genes <- c("AXL", "CD169", "SIGLEC1", "CD123", 
"HLA-DR", "CD45RA", "CD33", "CD5", "CX3CR1", "CD81", "CD1C", "CD141", "ITGAX",
"TCF4", "IL3RA")

signatures <- c("AMHD1", "BST2", "CD86", "CD80", "CX3CR1", "MX1", "RSAD2", 
"CXCL10", "CCL8", "CCR7", "SIGLEC1", "SIGLEC6", "AXL", "IL12A", "IL12B", "CD123", 
"IL6", "IL8", "CD274", "CD40", "CGAS", "TMEM173", "TLR3", "TLR7", "TLR9", "IFITM1", 
"IFITM2", "IFITM3")

HIV_macrophages_Decalf <- Decalf_genes <- c("OASL", "ISG15", "CXCL10", "CCL8", 
"HERC5", "IFI44L", "IFI44", "IFIT2", "CMPK2", "IFIT1", "IFIT3", "RSAD2", "CXCL11", 
"STAP1", "IFITM1", "TNFSF10", "GBP4")

pDC_signatures_basal <- c("PACSIN1", "IFIT2", "GZMB", 
"NOTCH4", "NLRP7", "IGF2R", "TLR7", "MYB", "IL28RA", "LAMP1", "XBP1", "CLEC4C", 
"RAB38", "TRAF3", "TCF4", "CD36", "NRP1", "TRAF4", "IRF7", "SEC61G", "NOTCH1", 
"LTB", "CD99", "CD68", "SEC61B", "STAT2", "FAM13A", "TNFRSF11B")
pDC_signatures_basal_See <- pDC_signatures_basal

preDC_signatures_basal <- c("ITGAX", "CD5", "CLEC10A", 
"RAB31", "ADAM8", "CD22", "CLEC12A", "CXCL16", "CD33", "IFITM3", "CEBPA", "CLEC4A", 
"SIGLEC6", "CD86", "THBD", "NLRP3", "CSF3R", "IL1R2", "CD48", "SEC14L1", "ITGA5", 
"RAB37", "IFITM2", "RUNX3", "KLF4", "LILRB3", "MYD88", "ICAM3", "ZBTB46", "FAM110A", 
"FAM38A", "CD93", "RAB24", "IL13RA1", "KLF2", "CD2", "SIGLEC10", "HLA-DOB", "CD244", 
"ICAM4", "FAM109A", "CD72", "BCR", "FAM89B", "TNFSF12", "RAB32", "GAPDH", "RAB8B", 
"CD44", "PACSIN2", "LY6G6C", "LTBR", "CD63", "KLF8", "HLA-DOA", "RAB7A", "FAM49B", 
"SEC11A", "HLA-DMB", "HLA-DRA", "FCGR3A")
preDC_signatures_basal_See <- preDC_signatures_basal

pDC_signatures_villani <- c("GZMB", "SELPLG", "DERL3", "PTPRCAP", "BCL11A", 
"LAMP5", "SLA2", "SELS", "NRP1", "SIDT1", "TCF4", "SLC15A4", "PTCRA", "IRF7",
"TRAF4")

meta <- c("nGene", "nUMI", "perc.mito", "perc.ribo")

tnf_alpha <- "^TNFA" # all gene symbols starting with "TNFA"
ifn_beta <- "^IFNB"
ifn_alpha <- "^IFNA"




cl_res <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)




GenerateSampleNamesExp1 <- function () {

	all_names <- c("H12-mock", "H24-mock", "H12-HIV", "H24-HIV", "H12-CpGA", "H24-CpGA", "H24-HIVAztNvp")
	return(all_names)
}


ExtractSampleNamesMockHIVtreatExp1 <- function () {

	all_names <- GenerateSampleNamesExp1()
	return(all_names[c(1:4),7])
}


ExtractSampleNamesMock24HIVExp1 <- function () {

	all_names <- GenerateSampleNamesExp1()
	return(all_names[c(2,3,4)])
}


GenerateSampleNamesExp2 <- function () {

	all_names <- c("MOCK-24h", "preDC-D0", "preDC-AD8-24h", "preDC-AD8-plus-Vpx-24h", "cDC2-D0", "cDC2-AD8-24h", "cDC2-AD8-plus-Vpx-24h", "P480-NRH1")
	return(all_names)
}


ExtractSampleNamespreDCExp2 <- function () {

	all_names <- GenerateSampleNamesExp2()
	return(all_names[c(2,3,4)])
}


ExtractSampleNamesMock24HIVExp2 <- function () {

	all_names <- GenerateSampleNamesExp2()
	return(all_names[c(1,3)])
}


ExtractSampleNamesMock24HIVtreatExp2 <- function () {

	all_names <- GenerateSampleNamesExp2()
	return(all_names[c(1,3,4)])
}


ExtractSampleNamesPreDCMockExp2 <- function () {

	all_names <- GenerateSampleNamesExp2()
	return(all_names[1:4])
}


ExtractSampleNamesT0patientHIVExp2 <- function () {

	all_names <- GenerateSampleNamesExp2()
	return(all_names[c(2,8)])
}


GenerateSampleNamesMVA <- function () {
	
	names <- c("uninf-rep1", "MVA-rep1", "uninf-rep2", "MVA-rep2")
	return(names)
}


GenerateSampleNamesColon <- function () {

	names <- c("DC-distal", "DC-proximal", "Macro-distal", "Macro-proximal")
	return(names)
}


GenerateSampleColors <- function () {

	exp1 <- GenerateSampleNamesExp1 ()
	exp2 <- GenerateSampleNamesExp2 ()
	patient2 <- "180817-P32-exvivoDC"
	exp_mva <- GenerateSampleNamesMVA ()
	exp_colon <- GenerateSampleNamesColon ()

	base_colors <- c(brewer.pal(10, "Paired"), "#000000", "#FFFFFF", brewer.pal(8, "Set2"), "#555555", brewer.pal(4, "Dark2"))

	l <- list()
	l[as.character(exp1)] <- base_colors[c(1,2,3,4,5,6,8)]
	l[as.character(exp2)] <- base_colors[c(16,11,10,17,13,14,15,12)]
	l[patient2] <- "orange"
	l[as.character(exp_mva)] <- base_colors[22:25]
	l[as.character(exp_colon)] <- base_colors[c(22,25,24,23)]

	return(l)
}


GeneratePreDCSampleColors <- function() {

	# 1,2   mock 12-24h
	# 3,4   HIV+vpx 12-24h
	# 5,6   CpGA 12-24h
	# 7,8   HIV+vpx+azt 12-24h
	# 9,10  HIV 12-24h
	# 11    uninfected 0h
	# 12    patient
	
	colors <- c(brewer.pal(10, "Paired"), "#000000", "#111111")
	# show_col(colors)

	return(colors)
}


ExtractColorsExp1 <- function() {

	colors <- GeneratePreDCSampleColors()
	return(colors[c(1:6,8)])
}


ExtractColorsMock24HIVtreatExp1 <- function() {

	colors <- GeneratePreDCSampleColors()
	return(colors[c(2,3,4,8)])
}


ExtractColorsMock24HIVExp1 <- function() {

	colors <- GeneratePreDCSampleColors()
	return(colors[c(2,3,4)])
}


ExtractColorspreDCExp2 <- function() {

	colors <- GeneratePreDCSampleColors()
	return(colors[c(11,10,4)])
}


ExtractColorsPreDCMockExp2 <- function() {

	colors <- GeneratePreDCSampleColors()
	return(colors[c(2,11,10,4)])
}


ExtractColorsMock24HIVExp2 <- function () {

	colors <- GeneratePreDCSampleColors()
	return(colors[c(2,10)])
}


ExtractColorsMock24HIVtreatExp2 <- function () {

	colors <- GeneratePreDCSampleColors()
	return(colors[c(2,10,8)])
}


ExtractColorsT0patientHIVExp2 <- function () {

	colors <- GeneratePreDCSampleColors()
	return(colors[c(11,12)])
}


LoadObject <- function(object_file)
{
	l <- load(object_file)
	object <- eval(parse(text=l))
	return(object)
}


# possibly add to an existing object
LoadMultiMatrix <- function(files.list, names.list, out.dir = "", min.cells=3, min.genes=200, object = NULL) {

	files_list <- as.character(as.array(read.table(files.list)[,1]))
	names_list <- as.character(as.array(read.table(names.list)[,1]))

	if (length(files_list) != length(names_list)) {
		cat("<files_list> and <names_list> must have the same length\n")
		return(1)
	}

	n <- length(files_list)

	cat("Load the datasets...")

	ncells <- 0

	# load the dataset
	data <- Read10X(data = files_list[1])

	if (n == 1 & is.null(object)) {

		cat(paste("Number of genes: ", nrow(data), "\n", sep=""))
		cat(paste("Number of cells: ", ncol(data), "\n", sep=""))
	
		object <- CreateSeuratObject(raw.data = data, project = names_list[1], min.cells = min.cells, min.genes = min.genes)

		cat(paste("Number of genes after filtering: ", nrow(object@data), "\n", sep=""))
		cat(paste("Number of cells after filtering: ", ncol(object@data), "\n", sep=""))

	} else {

		# Create the seurat object
		# Returns a Seurat object with the raw data stored in
		# object@raw.data. object@data, object@meta.data, object@ident, also
		# initialized.
		start <- 1
		if (is.null(object)) {
			object <- CreateSeuratObject(raw.data = data, project = names_list[1])
			start <- start + 1
		}

		ncells <- ncells + ncol(object@data)
		if (n > start) {
			for (i in start:(n-1)) {
				# load the dataset
				data <- Read10X(data = files_list[i])
				# Create the seurat object
				obj2 <- CreateSeuratObject(raw.data = data, project = names_list[i])
				ncells <- ncells + ncol(obj2@data)
				# merge data and filter genes
				# Keep all genes expressed in >= 3 cells (~0.1% of the data). Keep all cells with at least 200 detected genes
				# By doing this at the end, gene filtering is performed based on their expression on the cells of the combined datasets
				if (i == start) {
					object <- MergeSeurat(object1 = object, object2 = obj2, add.cell.id1 = object@ident, add.cell.id2 = obj2@ident, project="DC", do.normalize=FALSE)
				} else {
					object <- MergeSeurat(object1 = object, object2 = obj2, add.cell.id2 = obj2@ident, project="DC", do.normalize=FALSE)
				}
			}
		}
		# load the dataset
		data <- Read10X(data = files_list[n])
		# Create the seurat object
		obj2 <- CreateSeuratObject(raw.data = data, project = names_list[n])
		ncells <- ncells + ncol(obj2@data)

		cat("done\n")

		cat(paste("Number of genes: ", nrow(object@data), "\n", sep=""))
		cat(paste("Number of cells: ", ncells, "\n", sep=""))

		cat("Merge into a single object and filter (min.cells = ", min.cells, ", min.genes = ", min.genes, ")...\n")

		object <- MergeSeurat(object1 = object, object2 = obj2, add.cell.id1 = object@ident, add.cell.id2 = obj2@ident, project="DC", do.normalize=FALSE, min.cells = min.cells, min.genes = min.genes)

		cat("done\n")

		cat(paste("Number of genes after filtering: ", nrow(object@data), "\n", sep=""))
		cat(paste("Number of cells after filtering: ", ncol(object@data), "\n", sep=""))
	}

	return(object)
}


LoadMatrixSelectBarcodes <- function(matrix.dir, sample.id, out.dir = "", min.cells = 3, min.genes = 200, barcodes = NULL) {

	cat("Load the dataset...")

	# load the dataset
	data <- Read10X(data = matrix.dir)

	if (!is.null(barcodes))
	{
		# here select the input barcodes
		which_barcodes <- which(data@Dimnames[[2]] %in% barcodes)
		data <- data[, which_barcodes, drop = FALSE]
	}

	cat(paste("Number of genes: ", nrow(data), "\n", sep=""))
	cat(paste("Number of cells: ", ncol(data), "\n", sep=""))
	
	object <- CreateSeuratObject(raw.data = data, project = sample.id, min.cells = min.cells, min.genes = min.genes)

	cat(paste("Number of genes after filtering: ", nrow(object@data), "\n", sep=""))
	cat(paste("Number of cells after filtering: ", ncol(object@data), "\n", sep=""))

	return(object)
}


# Filter out the cells that exhibit a number of genes that either:
# 1. is lower than max(min.genes, local.min, mean-2*sd), where local.min is the estimated minimum between the two modes
# 2. is higher than min(max.genes, mean+2*sd)
# where mean and sd are computed on the distribution of the number of genes
BimodalnGeneFilter <- function(object, bimod = TRUE, outliers = FALSE, sd_factor = 2, perc.mito = 0, min.genes = 200, max.genes = MAXIMUM, 
						second_max_interval_start = 200, second_max_interval_end_subtract = 1) {

	local_min <- min.genes
	outlier <- max.genes

	cat(paste("Number of cells before nGene filtering: ", ncol(object@data), "\n", sep=""))

	if (bimod) {
		# Filter cells based on sample-specific nGene count
		d <- density(object@meta.data$nGene)
		first_max <- optimize(approxfun(d$x,-d$y), interval=c(min.genes, max.genes))$minimum
		second_max <- optimize(approxfun(d$x,-d$y), interval=c(second_max_interval_start, first_max-second_max_interval_end_subtract))$minimum
		# local minimum between the two modes
		local_min <- max(local_min, optimize(approxfun(d$x,d$y), interval=c(min(first_max,second_max), max(first_max,second_max)))$minimum)
	}

	if (outliers) {
		mean_nGene <- mean(object@meta.data$nGene)
		sd_nGene <- sd(object@meta.data$nGene)
		outlier <- mean_nGene + sd_factor*sd_nGene
#		local_min <- max(local_min, mean_nGene - sd_factor*sd_nGene)
	}

	object <- FilterCells(object, subset.names = c("nGene"), low.thresholds = c(local_min), high.thresholds = c(outlier))
	cat(paste("Number of cells after nGene filtering: ", ncol(object@data), "\n", sep=""))

	if (perc.mito > 0) {
		# Further filter cells on fraction of mito genes 
		object <- FilterCells(object, subset.names = c("perc.mito.UMI"), low.thresholds = c(-Inf), high.thresholds = c(perc.mito))
		cat(paste("Number of cells after filtering by mitochondrial gene expression: ", ncol(object@data), "\n", sep=""))
	}

	return(object)
}


# remove cells with high expression of either one of the genes in the input array
FilterCellsByHighGeneExpression <- function(object, genes, sd_factor = 3) {

	# compute the mean and sdev at the beginning
	high_thresholds <- low_thresholds <- rep(0,length(genes))
	for (i in 1:length(genes)) {
		if (genes[i] %in% rownames(object@data)) {
			m <- mean(object@data[genes[i],])
			s <- sd(object@data[genes[i],])
			high_thresholds[i] <- m + sd_factor*s	
		}
	}	

	for (i in 1:length(genes)) {
		cat("Filter cells based on expression of gene", genes[i], "\n")
		if (genes[i] %in% rownames(object@data)) {
			cat("Gene", genes[i], "is found in the list\n")
			object <- FilterCells(object, subset.names = genes[i], low.thresholds = low_thresholds[i] - EPSILON, high.thresholds = high_thresholds[i] + EPSILON)
			cat(paste("Number of cells after ", genes[i], " filtering: ", ncol(object@data), "\n", sep=""))
		}
	}
	
	return(object)
}


# remove cells with high expression of either one of the genes in the input array
FilterCellsByLowGeneExpression <- function(object, genes, sd_factor = 3) {

	high_thresholds <- low_thresholds <- rep(MAXIMUM, length(genes))
	for (i in 1:length(genes)) {
		if (genes[i] %in% rownames(object@data)) {
			m <- mean(object@data[genes[i],])
			s <- sd(object@data[genes[i],])
			low_thresholds[i] <- max(0, m - sd_factor*s)
		}
	}

	for (i in 1:length(genes)) {
		if (genes[i] %in% rownames(object@data)) {	
			object <- FilterCells(object, subset.names = genes[i], low.thresholds = low_thresholds[i] - EPSILON, high.thresholds = high_thresholds[i] + EPSILON)
			cat(paste("Number of cells after ", genes[i], " filtering: ", ncol(object@data), "\n", sep=""))
		}
	}
	
	return(object)
}


FilterCellsByGeneAuc <- function(object, genes, aucMaxRankPerc = 0.2, auc.cutoff = 0.25) {

	exprMatrix <- as.matrix(object@data)
	cells_rankings <- AUCell_buildRankings(exprMatrix)
	cells_AUC <- AUCell_calcAUC(genes, cells_rankings, aucMaxRank=nrow(cells_rankings)*aucMaxRankPerc)

	M <- as.matrix(getAUC(cells_AUC))
	for (i in 1:nrow(M)) {
		meta <- data.frame(M[i,])	
		colnames(meta) <- rownames(getAUC(cells_AUC))[i]
		rownames(meta) <- colnames(getAUC(cells_AUC))
		object <- AddMetaData(object = object, metadata = meta)
	}
	object <- FilterCells(object, subset.names = names(genes), low.thresholds = rep(0 - EPSILON, length(genes)), high.thresholds = rep(auc.cutoff, length(genes)))
	cat(paste("Number of cells after AUC filtering: ", ncol(object@data), "\n", sep=""))

	return(object)
}


PlotMetaData <- function (object, violin_plot_features, OUT_DIR, filtered = FALSE) {

	for (feature in violin_plot_features) {	
		if (!filtered) {
			plot.name <- paste(OUT_DIR, "/", feature, ".pdf", sep="")
		} else {
			plot.name <- paste(OUT_DIR, "/", feature, "_filtered.pdf", sep="")
		}
		pdf(plot.name, width=6, height=6)
		g <- VlnPlot(object = object, features.plot = feature)
		print(g)
		dev.off()
	}

	return()
}


NormScaleMatrix <- function(object, model="linear", vars.to.regress="nUMI") { 

	cat("Normalize and log-transform...")

	# LOG-NORMALIZE 
	# normalizes the gene expression measurements for each cell by the total expression, multiplies 
	# this by a scale factor (10,000 by default), and log-transforms the result
	obj <- NormalizeData(object = object, normalization.method = "LogNormalize", scale.factor = 10000)

	cat("done\n")

	cat("Scale and center the UMI counts...")

	# SCALING
	# regressing nUMI from obj@data
	# scale and center umi counts
	# saves the results in obj@scale.data 
	obj <- ScaleData(object = obj, vars.to.regress = vars.to.regress, model.use=model)

	cat("done\n")

	return(obj)
}


PCA <- function(object, out.dir, pcs.compute = 20, only_genes = NULL, suffix = "allgenes", do.jack = FALSE) {

	if (!file.exists(out.dir)){
		dir.create(file.path(getwd(), out.dir))
	}

	obj <- object
	reduction.name <- paste("pca", suffix, sep="_")

	cat("Run PCA...")

	# PCA
	sink(file = paste(out.dir, "/", reduction.name, ".txt", sep = ""))
	obj <- RunPCA(object = obj, pc.genes = only_genes, pcs.compute = pcs.compute)
	sink()

	cat("done\n")

	cat("Generate PCA dotplots...")

	ident <- as.vector(unique(object@ident)) # to assign the correct colors
	# NEW: assign the colors based on a general palette
	base_col <- GenerateSampleColors ()
	col <- unlist(base_col[match(ident,names(base_col))])

	# plot the components
	plot.title <- paste(out.dir, "/", reduction.name, ".pdf", sep="")
	pdf(plot.title, width=21, height=30)
	ncol <- 4
	nrow <- 7
	title <- reduction.name
	grid.newpage()
	pushViewport(viewport(layout = grid.layout(nrow = nrow+1, ncol = ncol, heights = unit(c(0.03, rep((1-0.05)/nrow, nrow)), "npc"))))
	grid.text(title, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:ncol, gp = gpar(fontsize = 35)))
	# stdev explained by the first i components, for i=1,...,pcs.compute
	stdev <- 0
	for (i in 1:floor(pcs.compute/2)) {
#		g <- DimPlot(object = obj, dim.1 = 2*i-1, dim.2 = 2*i, cols.use = colors[match(ident, all_names)], do.return = TRUE, reduction.use="pca", label.size=4)
		stdev1 <- stdev + obj@dr$pca@sdev[2*i-1]
		stdev2 <- stdev1 + obj@dr$pca@sdev[2*i]
		title <- paste("stdev(pc", 2*i-1, ") = ", round(stdev1, digits=1), ", stdev(pc", 2*i, ")=", round(stdev2, digits=1), "\n", sep="")
		g <- DimPlot(object = obj, dim.1 = 2*i-1, dim.2 = 2*i, cols.use = col, do.return = TRUE, label.size=4) + ggtitle(title)
		print(g, vp = viewport(layout.pos.row = ((i-1)%%(ncol*nrow))%/%ncol+2, layout.pos.col = (i-1)%%ncol+1))
		if (i%%(nrow*ncol) == 0 && floor(pcs.compute/2) > nrow*ncol) {
			grid.newpage()
			pushViewport(viewport(layout = grid.layout(nrow = nrow, ncol = ncol)))
		}
		stdev <- stdev2
	}
	dev.off()

	cat("done\n")

	### NB: The functions used below have reduction.type hardcoded ("pca")
	# So the correct pca ID is assigned at the end

	if (do.jack) {
		cat("Generate PCA JackStraw plot...")
	
		# Significant components
		obj <- JackStraw(object = obj, num.pc = pcs.compute, num.replicate = 100)
		plot.title <- paste(out.dir, "/pcaJack_", suffix, ".pdf", sep="")
		pdf(plot.title, width = 7, height = 7*ceiling(pcs.compute/25))
		print(JackStrawPlot(object = obj, nCol = 5, PCs = 1:pcs.compute))
		dev.off()

		cat("done\n")
	}

	cat("Generate PCA elbow plot...")

	# PCA stdev plot
	pdf(paste(out.dir, "/pcastdev_", suffix, ".pdf", sep=""))
	print(PCElbowPlot(object = obj, num.pc = pcs.compute))
	dev.off()

	cat("done\n")

	names(obj@dr)[names(obj@dr) == "pca"] <- reduction.name

	return(obj)
}


CellsClusters <- function(object, out.dir, dim = 1, res = cl_res, prefix = NULL, suffix = "allgenes", reduction.type = NULL, genes.use = NULL) {

	cat("Cluster cells from PCA...\n")

	obj <- object
	dims.use <- 1:dim

	# NB: FindClusters will use the "pca" reduction by default
	# this step ensures that "pca" is disabled whenever a set of genes is provided
	if (is.null(genes.use)) {
		if (is.null(reduction.type)) {
			cat("Either reduction.type or genes.use must be specified\n")
			return(obj)
		} else { 
			if (!any(names(obj@dr) == reduction.type)) {
				cat(paste("Reduction ", reduction.type, " has not been computed yet\n", sep=""))
				return(obj)
			}
		}
	} else { # !(is.null(genes.use)) 
		reduction.type <- NULL
		dims.use <- NULL
	}

	n <- 1
	for (r in res) {

		# Cells clustering
		cat(paste("(dim=", dim, ", resolution=", r, ")\n", sep=""))
		if (n == 1) {
			obj <- FindClusters(object = obj, reduction.type = reduction.type, dims.use = dims.use, resolution = r, print.output = 0, save.SNN = TRUE, force.recalc=TRUE, genes.use = genes.use)
		} else {
			# obj <- FindClusters(object = obj, reduction.type = reduction.type, dims.use = dims.use, resolution = r, print.output = 0, reuse.SNN = TRUE, genes.use = genes.use)
			obj <- FindClusters(object = obj, resolution = r, print.output = 0, reuse.SNN = TRUE, force.recalc=TRUE)
		}
		cat("done\n")

		if (!is.null(prefix)) {
			old_name <- paste("res.", r, sep="")
			new_name <- paste(prefix, ".", old_name, sep="")
			colnames(obj@meta.data)[colnames(obj@meta.data) == old_name] <- new_name
		}

#		write.table(obj@ident, file=paste(out.dir, "/cl_", suffix, "_dim", dim, "_res", r, ".tsv", sep=""), sep="\t")
		n <- n + 1
	}

	obj <- SetAllIdent(object = obj, id = "orig.ident")
	return(obj)
}


ClustersID <- function(res, prefix = NULL) {

	# id is the name of the cells cluster identity
	id <- paste("res.", res, sep="")
	if (!is.null(prefix)) {
		id <- paste(prefix, ".", id, sep="")
	}

	return(id)
}


SetClustersIdent <- function(object, res, prefix = NULL) {

	# id is the name of the cells cluster identity
	id <- ClustersID(res, prefix)
	if (!any(names(object@meta.data) == id)) {
		cat(paste("Clustering resolution ", res, " has not been computed yet [SetClustersIdent]\n", sep=""))
	} else {
		object <- SetAllIdent(object = object, id = id)
	}

	return(object)
}


# returns the number of cells where the (set of) gene(s) is detected
DetectedGenesPerCell <- function(object, id, genes) {

	s <- sum(object@raw.data[rownames(object@data) %in% genes,which(object@ident == id)] == 1)
	return(s)
}


# returns the number of cells where the (set of) gene(s) is detected
DetectedUMIsPerCell <- function(object, id, genes) {

	s <- sum(object@raw.data[rownames(object@data) %in% genes,which(object@ident == id)])
	return(s)
}


Silhouette <- function(object, data, out.prefix, res = 0.1*(1:8), prefix = NULL) {
	
	dist <- dist(data)
	for (r in res) {
		object <- SetClustersIdent(object, res = r, prefix)
		silh <- silhouette(as.numeric(as.matrix(object@ident)), dist)
		sink(paste(out.prefix, "_silhouette_", r, ".txt", sep=""))
		print(summary(silh))
		sink()
		pdf(paste(out.prefix, "_silhouette_", r, ".pdf", sep=""))
		plot(silh)
		dev.off()
	}

	return()
}

# remove the cells that belong to the specified clusters
# clusters are specified as a comma-separated string
# return the list of selected cells 
# select = TRUE will select the specified identities instead of removing them
FilterOutClusters <- function (object, id, cl_filter, select = FALSE) {

	if (cl_filter == "") {
		return(colnames(object@data))
	}

	cl <- unlist(strsplit(cl_filter, split=","))	

	if (select) {
		filtered_cells <- names(object@ident[object@meta.data[[id]] %in% cl])
	} else {
		filtered_cells <- names(object@ident[!(object@meta.data[[id]] %in% cl)])
	}

	return(filtered_cells)
}


ClusterNames <- function (object, id, exclude = NULL) {

	new.cl.names <- unique(sort(object@meta.data[[id]]))
	if (!is.null(exclude)) {
		new.cl.names <- new.cl.names[!(new.cl.names %in% exclude)]
	}

	return(new.cl.names)
}


# reassign cluster names
# renaming is specified as a comma-separated string and refers to clusters 0,1,2... in this order
# specify the array of clusters to be masked (i.e. exlcuded from the renaming
RenameClusters <- function (object, id, renaming = "", cl_filter = NULL) {

	old.cl.names <- orig.cl.names <- ClusterNames(object, id)
	if (!is.null(cl_filter)) {
		cl <- unlist(strsplit(cl_filter, split=","))
		old.cl.names <- orig.cl.names[!(orig.cl.names %in% cl)]
		object@meta.data[[id]] <- plyr::mapvalues(x = object@meta.data[[id]], from = cl, to = rep("X", length(cl))) # assign X to removed clusters
	}
	if (renaming != "") {
		new.cl.names <- unlist(strsplit(renaming, split=","))
		object@meta.data[[id]] <- plyr::mapvalues(x = object@meta.data[[id]], from = old.cl.names, to = new.cl.names) # rename only filtered clusters according to the mapping
	}

	return(object)
}


RuntSNEandPlot <- function (object, dim, seed, reduct = "tsne", pdf.plot, pt.size = 1, label.size = 12) {

	obj <- RunTSNE(object = object, dims.use = 1:dim, do.fast = TRUE, seed.use = seed, reduction.name = reduct)
	title <- paste("dim = ", dim, ", seed = ", seed, sep = "")	

	pdf(pdf.plot, width=5, height=5)
	DimPlot(object = obj, no.legend=TRUE, do.label=TRUE, label.size=label.size, no.axes=TRUE, pt.size = pt.size, plot.title = title, reduction.use = reduct)
	dev.off()
}


# names.list is a reordering of ident (for legend only)
# tSNEplotSamples <- function(object, out.dir, title, names.list, seed = 1, pt.size = 1, reduction.use = "tsne", cells.use = NULL) {
tSNEplotSamples <- function(object, out.dir, title, seed = 1, pt.size = 1, reduction.use = "tsne", cells.use = NULL) {

	# tSNE plot

#	names <- as.character(as.array(read.table(names.list)[,1]))

#	ident <- as.vector(unique(object@ident)) # to assign the correct colors
	ident <- levels(object@ident)

	# NEW: assign the colors based on a general palette
	base_col <- GenerateSampleColors ()
	col <- unlist(base_col[match(ident,names(base_col))])

	cat("Generate tSNE plot...")

	pdf(paste(out.dir, "/", reduction.use , "_seed", seed, "_samples.pdf", sep=""), width=5, height=5)
#	DimPlot(object = object, no.legend = TRUE, no.axes = TRUE, pt.size = pt.size, plot.title = title, cols.use = colors[match(ident, all_names)], reduction.use = reduction.use, cells.use = cells.use)
	DimPlot(object = object, no.legend = TRUE, no.axes = TRUE, pt.size = pt.size, plot.title = title, cols.use = col, reduction.use = reduction.use, cells.use = cells.use)
	dev.off()

	if (length(ident) > 1) {
		pdf(paste(out.dir, "/", reduction.use, "_seed", seed, "_samples_legend.pdf", sep = ""))
		plot(x = 0, y = 0, pch = NA, axes = FALSE, ann = FALSE)
	#	legend(0, 1, col = colors[match(names, all_names)], legend = names, bty = "n", lwd = 20, seg.len = 0.1, y.intersp = 1.5, x.intersp = 1.5)
		legend(0, 1, col = col, legend = ident, bty = "n", lwd = 20, seg.len = 0.1, y.intersp = 1.5, x.intersp = 1.5)
		dev.off()
	}

	return()
}


# names.list is a reordering of ident (for legend only)
# Create a multiple plot containing the tSNE for each seed and each perplexity
tSNEplotSamplesMulti <- function(object, out.dir, seed, perplexity, pt.size = 1, reduction.use = "tsne", cells.use = NULL) {

	# tSNE plot
	ident <- as.vector(unique(object@ident)) # to assign the correct colors

	# NEW: assign the colors based on a general palette
	base_col <- GenerateSampleColors ()
	col <- unlist(base_col[match(ident,names(base_col))])

	cat("Generate tSNE plot...")

	plot.name = paste(out.dir, "/", reduction.use , "_samples.pdf", sep="")
	pdf(plot.name, width=21, height=30)
	ncol <- 3
	nrow <- 5
	title <- reduction.use
	grid.newpage()
	pushViewport(viewport(layout = grid.layout(nrow = nrow+1, ncol = ncol, heights = unit(c(0.03, rep((1-0.05)/nrow, nrow)), "npc"))))
	grid.text(title, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:ncol, gp = gpar(fontsize = 35)))
	n <- 1
	for (i in 1:length(seed)) {
		s <- seed[i]
		for (j in 1:length(perplexity)) {
			p <- perplexity[j]
			red <- paste(paste(reduction, ".seed", s, ".perp", p, sep=""))
			plot_title <- paste(paste("seed = ", s, ", perp. = ", p, sep=""))
			g <- DimPlot(object = object, no.axes = TRUE, pt.size = pt.size, cols.use = col, reduction.use = red, cells.use = cells.use, do.return = TRUE) + ggtitle(plot_title)
			print(g, vp = viewport(layout.pos.row = ((n-1)%%(ncol*nrow))%/%ncol+2, layout.pos.col = (n-1)%%ncol+1))
			if (n%%(nrow*ncol) == 0 && n > nrow*ncol) {
				grid.newpage()
				pushViewport(viewport(layout = grid.layout(nrow = nrow+1, ncol = ncol, heights = unit(c(0.03, rep((1-0.05)/nrow, nrow)), "npc"))))
				grid.text(title, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:ncol, gp = gpar(fontsize = 35)))
			}
			n <- n + 1
		}
	}
	dev.off()

	return()
}


# names.list is a reordering of ident (for legend only)
# Create a multiple plot containing the tSNE for each seed and each perplexity
UMAPplotSamplesMulti <- function(object, out.dir, seed, neighbors, min_dist, pt.size = 1, reduction.use = "tsne", cells.use = NULL) {

	# tSNE plot
	ident <- as.vector(unique(object@ident)) # to assign the correct colors

	# NEW: assign the colors based on a general palette
	base_col <- GenerateSampleColors ()
	col <- unlist(base_col[match(ident,names(base_col))])

	cat("Generate UMAP plot...")

	plot.name = paste(out.dir, "/", reduction.use , "_samples.pdf", sep="")
	pdf(plot.name, width=21, height=30)
	ncol <- 3
	nrow <- 5
	title <- reduction.use
	grid.newpage()
	pushViewport(viewport(layout = grid.layout(nrow = nrow+1, ncol = ncol, heights = unit(c(0.03, rep((1-0.05)/nrow, nrow)), "npc"))))
	grid.text(title, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:ncol, gp = gpar(fontsize = 35)))
	n <- 1
	for (i in 1:length(seed)) {
		s <- seed[i]
		for (j in 1:length(neighbors)) {
			nn <- neighbors[j]
			for (k in 1:length(min_dist)) {
				d <- min_dist[k]
				red <- paste(paste(reduction, ".seed", s, ".neigh", nn, ".dist", d, sep=""))
				plot_title <- paste(paste("seed = ", s, ", neigh. = ", nn, ", min. dist. = ", d, sep=""))
				g <- DimPlot(object = object, no.axes = TRUE, pt.size = pt.size, cols.use = col, reduction.use = red, cells.use = cells.use, do.return = TRUE) + ggtitle(plot_title)
				print(g, vp = viewport(layout.pos.row = ((n-1)%%(ncol*nrow))%/%ncol+2, layout.pos.col = (n-1)%%ncol+1))
				if (n%%(nrow*ncol) == 0 && n > nrow*ncol) {
					grid.newpage()
					pushViewport(viewport(layout = grid.layout(nrow = nrow+1, ncol = ncol, heights = unit(c(0.03, rep((1-0.05)/nrow, nrow)), "npc"))))
					grid.text(title, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:ncol, gp = gpar(fontsize = 35)))
				}
				n <- n + 1		
			}
		}
	}
	dev.off()

	return()
}

# id = id assigned to each cell (e.g. after clustering)
# id <- paste("cl_res", r, sep="")
tSNEplotClusters <- function(object, out.dir, title, seed = 1, res, prefix = NULL, new.cluster.ids = NULL, pt.size = 1, reduction.use = "tsne", cells.use = NULL, label.size = 12) {

	for (r in res) {
		
		object <- SetClustersIdent(object, r, prefix)

		if (length(new.cluster.ids) > 0) {
			current.cluster.ids <- unique(sort(object@ident))
			object@ident <- plyr::mapvalues(x = object@ident, from = current.cluster.ids, to = new.cluster.ids)
		}
	
		plot.title <- paste(title, " res=", r, sep = "")

		pdf(paste(out.dir, "/", reduction.use , "_cl_res", r, ".pdf", sep=""), width=5, height=5)
		g <- DimPlot(object = object, no.legend=TRUE, do.label=TRUE, no.axes=TRUE, pt.size = pt.size, reduction.use = reduction.use, cells.use = cells.use, do.return = TRUE, label.size = label.size)
		g <- g + ggtitle(plot.title)
		print(g)
		dev.off()
	}
}


# names.list is a reordering of ident (for legend only)
# Create a multiple plot containing the tSNE for each seed and each perplexity
tSNEplotClustersMulti <- function(object, out.dir, seed, perplexity, res, prefix = NULL, pt.size = 1, reduction.use = "tsne", cells.use = NULL, label.size = 8) {

	cat("Generate tSNE plot...")

	plot.name = paste(out.dir, "/", reduction.use , "_clusters.pdf", sep="")
	pdf(plot.name, width=21, height=30)
	ncol <- 4
	nrow <- 5
	title <- reduction.use
	grid.newpage()
	pushViewport(viewport(layout = grid.layout(nrow = nrow+1, ncol = ncol, heights = unit(c(0.03, rep((1-0.05)/nrow, nrow)), "npc"))))
	grid.text(title, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:ncol, gp = gpar(fontsize = 35)))
	n <- 1
	for (k in 1:length(res)) {
		
		object <- SetClustersIdent(object, res[k], prefix)

		for (i in 1:length(seed)) {
			s <- seed[i]
			for (j in 1:length(perplexity)) {
				p <- perplexity[j]
				red <- paste(paste(reduction, ".seed", s, ".perp", p, sep=""))
				plot_title <- paste(paste("seed = ", s, ", perp. = ", p, ", res = ", res[k], sep=""))
				g <- DimPlot(object = object, no.legend=TRUE, do.label=TRUE, no.axes = TRUE, label.size=label.size, pt.size = pt.size, reduction.use = red, cells.use = cells.use, do.return = TRUE) + ggtitle(plot_title)
				print(g, vp = viewport(layout.pos.row = ((n-1)%%(ncol*nrow))%/%ncol+2, layout.pos.col = (n-1)%%ncol+1))
				if (n%%(nrow*ncol) == 0 && n < length(res)*length(seed)*length(perplexity)) {
					grid.newpage()
					pushViewport(viewport(layout = grid.layout(nrow = nrow+1, ncol = ncol, heights = unit(c(0.03, rep((1-0.05)/nrow, nrow)), "npc"))))
					grid.text(title, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:ncol, gp = gpar(fontsize = 35)))
				}
				n <- n + 1
			}
		}
	}

	dev.off()

	return()
}


# names.list is a reordering of ident (for legend only)
# Create a multiple plot containing the UMAP for each seed and each perplexity
UMAPplotClustersMulti <- function(object, out.dir, seed, neighbors, min_dist, res, prefix = NULL, pt.size = 1, reduction.use = "tsne", cells.use = NULL, label.size = 8) {

	cat("Generate UMAP plot...")

	plot.name = paste(out.dir, "/", reduction.use , "_clusters.pdf", sep="")
	pdf(plot.name, width=21, height=30)
	ncol <- 4
	nrow <- 5
	title <- reduction.use
	grid.newpage()
	pushViewport(viewport(layout = grid.layout(nrow = nrow+1, ncol = ncol, heights = unit(c(0.03, rep((1-0.05)/nrow, nrow)), "npc"))))
	grid.text(title, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:ncol, gp = gpar(fontsize = 35)))
	n <- 1
	for (k in 1:length(res)) {
		
		object <- SetClustersIdent(object, res[k], prefix)

		for (i in 1:length(seed)) {
			s <- seed[i]
			for (j in 1:length(neighbors)) {
				nn <- neighbors[j]
				for (h in 1:length(min_dist)) {
					d <- min_dist[h]
					red <- paste(paste(reduction, ".seed", s, ".neigh", nn, ".dist", d, sep=""))
					plot_title <- paste(paste("seed = ", s, ", neigh. = ", nn, ", min.dist. = ", d, ", res = ", res[k], sep=""))
					g <- DimPlot(object = object, no.legend=TRUE, do.label=TRUE, no.axes = TRUE, label.size=label.size, pt.size = pt.size, reduction.use = red, cells.use = cells.use, do.return = TRUE) + ggtitle(plot_title)
					print(g, vp = viewport(layout.pos.row = ((n-1)%%(ncol*nrow))%/%ncol+2, layout.pos.col = (n-1)%%ncol+1))
					if (n%%(nrow*ncol) == 0 && n < length(res)*length(seed)*length(neighbors)*length(min_dist)) {
						grid.newpage()
						pushViewport(viewport(layout = grid.layout(nrow = nrow+1, ncol = ncol, heights = unit(c(0.03, rep((1-0.05)/nrow, nrow)), "npc"))))
						grid.text(title, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:ncol, gp = gpar(fontsize = 35)))
					}
					n <- n + 1
				}
			}
		}
	}

	dev.off()

	return()
}


# show: array of clusters IDs to show in the plot
# new.cluster.ids is neeeded either to rename or to select a subset
tSNEplotClusterSel <- function(object, plot.name, title, seed = 1, res, prefix = NULL, new.cluster.ids = NULL, pt.size = 1, reduction.use = "tsne", cells.use = NULL, cl.show = NULL, label.size = 8) {

	object <- SetClustersIdent(object, res, prefix)

	cl_show <- unlist(strsplit(cl.show, split=","))

	# select colors
	cluster.ids <- unique(sort(object@ident))
	palette <- hue_pal()(length(cluster.ids[!(cluster.ids %in% c("X"))]))
	cols.use <- rep("gray", length(cluster.ids))
	cols.use[cluster.ids %in% cl_show] <- palette[cluster.ids %in% cl_show]
	
	# reassign names (delete all but shown clusters' labels)
	names <- rep("", length(cluster.ids))
	names[cluster.ids %in% cl_show] <- cluster.ids[cluster.ids %in% cl_show]
	object@ident <- plyr::mapvalues(x = object@ident, from = cluster.ids, to = names)

	pdf(plot.name, width=5, height=5)
	DimPlot(object = object, no.legend=TRUE, do.label=TRUE, label.size=label.size, no.axes=TRUE, pt.size = pt.size, plot.title = title, reduction.use = reduction.use, cells.use = cells.use, cols.use = cols.use)
	dev.off()

	return()
}


GeneSynonymToSymbol <- function (genes, org="Hs") {

	genes_official <- alias2SymbolTable(genes, species=org) # to convert gene aliases to official gene names
	genes_official <- genes_official[!is.na(genes_official)]

	cat(paste(length(genes), "genes read,", length(genes_official), "gene official names\n"))

	return(genes_official)
}


# compute the fraction of transcripts of the input gene group and add this info to a metadata slot
# group_name identifies the gene group in the metadata slot
GeneGroupExpression <- function(object, genes, group_name, convert.to.official.name = FALSE, org="Hs") {
	
	gene_names <- rownames(object@raw.data)
	if (convert.to.official.name) {
		gene_names <- GeneSynonymToSymbol(gene_names, org=org)
		genes <- GeneSynonymToSymbol(genes, org=org)
	}

	# FIXME: this does not work if there is just one gene is expressed in some cells!!!
	sum_genes_UMI <- apply(object@raw.data[gene_names %in% genes,], 2, sum)
	sum_UMI <- apply(object@raw.data, 2, sum) # N.B. issue with large datasets here. FIXME: use sparse matrix implementation
	perc_genes_UMI <- sum_genes_UMI/sum_UMI
	log_perc_genes_UMI <- log10(10000*perc_genes_UMI+1)
	
	sum_nGene <- apply(object@raw.data[gene_names %in% genes,], 2, function(x) sum(x > 0))

	object <- AddMetaData(object = object, metadata = sum_genes_UMI, col.name = paste(group_name, "UMI", sep="."))
	object <- AddMetaData(object = object, metadata = perc_genes_UMI, col.name = paste("perc", group_name, "UMI", sep="."))
	object <- AddMetaData(object = object, metadata = log_perc_genes_UMI, col.name = paste("log.perc", group_name, "UMI", sep="."))
	object <- AddMetaData(object = object, metadata = sum_nGene, col.name = paste(group_name, "genes", sep="."))

	return(object)
}


HIVexpression <- function(object) {

	object <- GeneGroupExpression(object, HIV_genes, "HIV1")
	return(object)
}

# compute the fraction of HIV transcripts and add this info to a metadata slot
# HIVexpression <- function(object) {
#
#	# FIXME: this does not work if there is just one HIV gene expressed in some cells!!!
#	sum_HIV1_UMI <- apply(object@raw.data[rownames(object@raw.data) %in% HIV_genes,], 2, sum)
#	sum_UMI <- apply(object@raw.data, 2, sum) # N.B. issue with large datasets here. FIXME: use sparse matrix implementation
#	perc_HIV1_UMI <- sum_HIV1_UMI/sum_UMI
#	log_perc_HIV1_UMI <- log10(10000*perc_HIV1_UMI+1)
#
#	object <- AddMetaData(object = object, metadata = sum_HIV1_UMI, col.name = "HIV1.UMI")
#	object <- AddMetaData(object = object, metadata = perc_HIV1_UMI, col.name = "perc.HIV1.UMI")
#	object <- AddMetaData(object = object, metadata = log_perc_HIV1_UMI, col.name = "log.perc.HIV1.UMI")
#
#	return(object)
# }


# compute the fraction of HIV transcripts and add this info to a metadata slot
HIVexpressionBulk <- function(data) {

	sum_UMI <- colSums(data, sparseResult=TRUE)
	data <- data[rownames(data) %in% HIV_genes, , drop=FALSE]
	sum_HIV1_UMI <- colSums(data, sparseResult=TRUE)
	perc.HIV1.UMI <- sum_HIV1_UMI/sum_UMI

	return(perc.HIV1.UMI)
}


# compute the correlation between the expression levels of two arrays of gene names
# merge or not the expression values of the genes in the array
# return a matrix with correlation values
ExpressionCorrelation <- function(object, genes1, genes2, merge1 = FALSE, merge2 = FALSE) {

	nrow <- length(genes1)
	ncol <- length(genes2)
	rownames <- genes1
	colnames <- genes2
	if (merge1) {
		nrow <- 1
		rownames <- "set1"
		L <- list(as.character(as.matrix(genes1)))
		object <- AddModuleScore(object, genes.list = L)
		expr1 <- object@meta.data$Cluster1
	}
	if (merge2) {
		ncol <- 1
		colnames <- "set2"
		L <- list(as.character(as.matrix(genes2)))
		object <- AddModuleScore(object, genes.list = L)
		expr2 <- object@meta.data$Cluster1
	}
	Cor <- matrix(0, nrow = nrow, ncol = ncol)	
	rownames(Cor) <- rownames
	colnames(Cor) <- colnames
	i <- 1
	if (merge1) {
		j <- 1
		if (merge2) {
			Cor[i,j] <- cor(expr1, expr2)	
		} else {
			for (gene2 in genes2) {
				expr2 <- object@data[gene2,]
				Cor[i,j] <- cor(expr1, expr2)
				j <- j+1
			}
		}
	} else {
		for (gene1 in genes1) {
			expr1 <- object@data[gene1,]
			j <- 1
			if (merge2) {
				Cor[i,j] <- cor(expr1, expr2)	
			} else {
				for (gene2 in genes2) {
					expr2 <- object@data[gene2,]
					Cor[i,j] <- cor(expr1, expr2)	
					j <- j+1
				}
			}
			i <- i+1
		}
	}
	return(Cor)
}

# tSNE must be computed already
tSNEplotHIV <- function(object, out.dir, reduction.use = "tsne", seed = 1, cells.use = NULL, pt.size = 1, no.legend = TRUE) {

	# compute the fraction of HIV transcripts
	object <- HIVexpression (object)

	if (no.legend) {
		plot_name <- paste(out.dir, "/", reduction.use, "_HIV1.pdf", sep="")
	} else {
		plot_name <- paste(out.dir, "/", reduction.use, "_HIV1_withlegend.pdf", sep="")
	}
	pdf(plot_name, width=5+(!no.legend), height=5)
	FeaturePlot(object = object, features.plot = c("log.perc.HIV1.UMI") , cols.use = c("gray","blue"), 
		reduction.use = reduction.use, pt.size = pt.size, no.axes=TRUE, cells.use = cells.use, no.legend = no.legend)
	dev.off()

	return(object)
}


tSNEplotALL <- function(object, out.dir, names.list = "", names = "", title, seed = 1, res = cl_res, prefix = NULL, new.cluster.ids = NULL, pt.size = 1, reduction.use = "tsne", cells.use = NULL, label.size = 8) {

#	if (names.list == "" && names == "") {
#		names <- as.vector(unique(object@ident))
#	} else {
#		if (names == "") {
#			names <- as.character(as.array(read.table(names.list)[,1]))
#		}
#	}

	# tSNE plot coloured by sample
#	tSNEplotSamples(object = object, out.dir = out.dir, names = names, title = title, seed = seed, pt.size = pt.size, reduction.use = reduction.use, cells.use = cells.use) 
	tSNEplotSamples(object = object, out.dir = out.dir, title = title, seed = seed, pt.size = pt.size, reduction.use = reduction.use, cells.use = cells.use) 

	# tSNE plot coloured by cluster
	tSNEplotClusters(object = object, out.dir = out.dir, title = title, seed = seed, res = res, prefix = prefix, pt.size = pt.size, reduction.use = reduction.use, cells.use = cells.use, label.size = label.size)

	# tSNE plot coloured by the fraction of HIV-1 transcripts
	tSNEplotHIV(object = object, out.dir = out.dir, reduction.use = reduction.use, seed = seed, cells.use = cells.use, pt.size = pt.size)

	# tSNE plot coloured by the fraction of HIV-1 transcripts (with legend)
	tSNEplotHIV(object = object, out.dir = out.dir, reduction.use = reduction.use, seed = seed, cells.use = cells.use, pt.size = pt.size, no.legend = FALSE)

	return()
}


SamplesComposition <- function(object, out.dir, res, prefix = NULL, sample.names = NULL, cl.names = NULL) {

	orig_ident <- object@ident
	object <- SetClustersIdent(object, res, prefix)

	if (is.null(sample.names)) {
		sample.names <- unique(orig_ident)
	}

	if (is.null(cl.names)) {
		cl.names <- unique(object@ident)
	}

	cl_distribution <- matrix(NA, nrow = length(sample.names)*length(cl.names), ncol = 3)
	cl_distribution[,3] <- 0 

	n <- 1
	for (s in sample.names[length(sample.names)-0:(length(sample.names)-1)]) {
		for (c in cl.names) {
			count <- sum(orig_ident == s & object@ident == c)
			cl_distribution[n,] <- c(s, c, count)
			n <- n + 1
		}
	}

	cl_distribution <- as.data.frame(cl_distribution)
	colnames(cl_distribution) <- c("sample", "cluster", "num") 
	write.table(cl_distribution, paste(out.dir, "/samples_clusters_composition.", res, ".tsv", sep=""), sep="\t", row.names=FALSE, quote=FALSE)

	return()
}


BarplotSampleComposition <- function(out.dir, ordering = NULL) {

	# plot cluster distribution
	data <- read.table(paste(out.dir, "/samples_clusters_composition.", res, ".tsv", sep=""), sep="\t", header=TRUE)
	data$cluster <- as.factor(data$cluster)
	
	if (length(ordering) > 0) {
		data$sample <- factor(data$sample, levels = ordering)
	}

	pdf(paste(out.dir, "/samples_clusters_composition.", res, ".pdf", sep=""), width=6, height=3)
	print(ggplot(data, aes(x=sample, y=num, fill=cluster)) + geom_bar(stat="identity") + theme_minimal() + coord_flip() + xlab("") + ylab("Number of cells") + theme(axis.text.y=element_text(size=15)))
	dev.off()

	pdf(paste(out.dir, "/samples_clusters_composition_norm.", res, ".pdf", sep=""), width=6, height=3)
	print(ggplot(data, aes(x=sample, y=num, fill=cluster)) + geom_bar(stat="identity", position="fill") + theme_minimal() + coord_flip() + xlab("") + ylab("Fraction of cells") + theme(axis.text.y=element_text(size=15)))
	dev.off()
	
	return()
}


# Create a multiple plot containing the tSNE for each seed and each perplexity
BarplotSampleCompositionMulti <- function(object, out.dir, res, prefix = NULL, ordering = NULL) {

	cat("Generate sample composition plot...")

	plot.name = paste(out.dir, "/samples_clusters_composition.pdf", sep="")
	pdf(plot.name, width=21, height=30)
	ncol <- 3
	nrow <- 7
	title <- "Samples composition"
	grid.newpage()
	pushViewport(viewport(layout = grid.layout(nrow = nrow+1, ncol = ncol, heights = unit(c(0.03, rep((1-0.05)/nrow, nrow)), "npc"))))
	grid.text(title, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:ncol, gp = gpar(fontsize = 35)))
	n <- 1
	for (k in 1:length(res)) {
		
		# plot cluster distribution
		data <- read.table(paste(out.dir, "/samples_clusters_composition.", res[k], ".tsv", sep=""), sep="\t", header=TRUE)
		data$cluster <- as.factor(data$cluster)
		data$sample <- as.factor(data$sample)

		if (length(ordering) > 0) {
			data$sample <- factor(data$sample, levels = ordering)
		}
	
		if (n%%(nrow*ncol) == 1) {
			grid.text(title, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:ncol, gp = gpar(fontsize = 35)))
		}

		plot_title <- paste(paste("res = ", res[k], sep=""))
		g <- ggplot(data, aes(x=sample, y=num, fill=cluster)) + geom_bar(stat="identity") + theme_minimal() + coord_flip() 
		g <- g + xlab("") + ylab("Number of cells") + theme(axis.text.y=element_text(size=15)) + ggtitle(plot_title)
		print(g, vp = viewport(layout.pos.row = ((n-1)%%(ncol*nrow))%/%ncol+2, layout.pos.col = (n-1)%%ncol+1))

		n <- n + 1
	}

	dev.off()

	return()
}


# NB: this works just for human!!!
ComputeMeta <- function(object = object) {

	# mitochondrial genes
	mito <- grep(pattern = "^MT-", x = rownames(object@data), value = TRUE)
	percent.mito <- Matrix::colSums(object@raw.data[mito,])/Matrix::colSums(object@raw.data)
	object <- AddMetaData(object = object, metadata = percent.mito, col.name = "perc.mito")

	# ribosomal genes
	ribo <- grep(pattern = "(^RPL)|(^MRPL)|(^MRPS)|(^RPS)", x = rownames(object@data), value = TRUE)
	percent.ribo <- Matrix::colSums(object@raw.data[ribo,])/Matrix::colSums(object@raw.data)
	object <- AddMetaData(object = object, metadata = percent.ribo, col.name = "perc.ribo")

	return(object)
}


# generate a violin plot for each input gene in the list, grouped by cluster
# if res=0, keep the original identity instead (i.e. samples)
ViolinPlotClusters <- function(object, plot.name, markers, res = 0, prefix = NULL, clusters = NULL, title = "", check = TRUE) {

	# remove the signatures that are not found in the list of genes
	if (check) {
		markers <- markers[which(markers %in% object@data@Dimnames[[1]])]
	}

	if (res > 0) {
		object <- SetClustersIdent(object, res, prefix)
		if (length(clusters) == 0) {
			clusters <- levels(object@ident)
		}
	}

	nrow <- 4
	ncol <- 3

	if (title == "") {
		title <- "Gene expression"
	}
	
	pdf(plot.name, width=21, height=30)
	grid.newpage()
	pushViewport(viewport(layout = grid.layout(nrow = nrow+1, ncol = ncol, heights = unit(c(0.03, rep((1-0.05)/nrow,nrow)), "npc"))))
	grid.text(title, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:ncol, gp = gpar(fontsize = 35)))
	for (i in 1:length(markers)) {
		g <- VlnPlot(object = object, ident.include = clusters, features.plot = markers[i])
		print(g, vp = viewport(layout.pos.row = ((i-1)%%(ncol*nrow))%/%ncol+2, layout.pos.col = (i-1)%%ncol+1))
		if (i%%(nrow*ncol) == 0) {
			grid.newpage()
#			pushViewport(viewport(layout = grid.layout(nrow = nrow, ncol = ncol)))
			pushViewport(viewport(layout = grid.layout(nrow = nrow+1, ncol = ncol, heights = unit(c(0.03, rep((1-0.05)/nrow,nrow)), "npc"))))
			grid.text(title, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:ncol, gp = gpar(fontsize = 35)))
		}
	}
	dev.off()

	return()
}


# generate a violin plot for each input gene in the list, grouped by cluster
ViolinPlotClustersSingleMarker <- function(object, plot.name, marker, res, prefix = NULL, clusters = NULL, title = "", title.size = 30, check = TRUE) {

	# remove the signatures that are not found in the list of genes
	if (check) {
		if (!(marker %in% object@data@Dimnames[[1]])) {
			cat("Marker", marker, "not found among the genes\n", sep=" ")
			return()
		}
	}

	object <- SetClustersIdent(object, res, prefix)
	
	if (length(clusters) == 0) {
		clusters <- levels(object@ident)
	}

	if (title == "") {
		title <- marker
	}

	pdf(plot.name, width=6, height=6)
	g <- VlnPlot(object = object, ident.include = clusters, features.plot = marker) 
	g <- g + ggtitle(label=title) + theme(plot.title=element_text(size=title.size))
	print(g)
	dev.off()

	return()
}


ViolinPlotClustersClass <- function(object, plot.name, marker_pattern, res, prefix = NULL, clusters = NULL, title = "") {

	v <- rownames(object@raw.data)
	markers <- v[grep(marker_pattern, v)]

	ViolinPlotClusters(object, plot.name, markers, res, prefix = prefix, clusters = clusters, title = title)

	return()
}


# N.B.: here the clustering is performed *without* considering HIV genes
# Nevertheless, HIV genes might be found among the significant DEG
# NEW: compute also the average expression
# the logFC is defined as log(expr1+1) - log(expr2+1), where log is the natural logarithm
GeneMarkersTable <- function (object, out.name, ident.1, ident.2, test.use = "wilcox", min.pct = 0.1, logFC = 0.25) {

	cells.1 <- WhichCells(object = object, ident = ident.1)
	exp1 <- apply(object@data[, cells.1, drop = F], 1, function(x) mean(x = expm1(x = x)))

	cells.2 <- WhichCells(object = object, ident = ident.2)
	exp2 <- apply(object@data[, cells.2, drop = F], 1, function(x) mean(x = expm1(x = x)))

	cluster_markers <- FindMarkers(object = object, ident.1 = ident.1, ident.2 = ident.2, min.pct = min.pct, logfc.threshold = logFC, test.use = test.use)
		
	df <- data.frame(cluster_markers, avg_exp.1 = exp1[rownames(cluster_markers)], avg_exp.2 = exp2[rownames(cluster_markers)])
	write.table(df, file=out.name, sep="\t", quote=FALSE)

	return()
}


# NB: the row data are used automatically in FindMarkers() when a UMI counts-based method is selected!!!
ClusterGeneMarkersByPairs <- function(object, out.dir, res, prefix = NULL, clusters = NULL, test.use = "wilcox", min.pct = 0.1, logFC = 0.25) {

	object <- SetClustersIdent(object, res, prefix)

	if (length(clusters) == 0) {
#		clusters <- sort(unique(object@ident)) # all clusters # FIXED!! Use levels otherwise clusters are considered as factors
		clusters <- levels(object@ident)
	}

	for (i in 1:(length(clusters)-1)) {	
		for (j in (i+1):length(clusters)) {
			deg.table <- paste(out.dir, "/DEG_", test.use, "_cl", clusters[i], "-", clusters[j], "_res", res, ".tsv", sep="")
			GeneMarkersTable(object, out.name = deg.table, ident.1 = clusters[i], ident.2 = clusters[j], test.use = test.use, min.pct = min.pct, logFC = logFC)			
		}
	}

	return()
}


# NB: the row data are used automatically in FindMarkers() when a UMI counts-based method is selected!!!
ClusterGeneMarkersVsAll <- function(object, out.dir, res, prefix = NULL, clusters = NULL, test.use = "wilcox", min.pct = 0.1, logFC = 0.25) {

	object <- SetClustersIdent(object, res, prefix)

	if (length(clusters) == 0) {
#		clusters <- sort(unique(object@ident)) # all clusters
		clusters <- levels(object@ident)
	}

	for (i in 1:length(clusters)) {
		deg.table <- paste(out.dir, "/DEG_", test.use, "_cl", clusters[i], "-all_res", res, ".tsv", sep="")
		GeneMarkersTable(object, out.name = deg.table, ident.1 = clusters[i], ident.2 = setdiff(clusters, c(clusters[i])), test.use = test.use, min.pct = min.pct, logFC = logFC)	
	}

	return()
}


# NEW: compute also the average expression
# the logFC is defined as log(expr1+1) - log(expr2+1), where log is the natural logarithm
# NEW TABLE FORMAT
GeneMarkersTableNEW <- function(object, out.name, ident.1, ident.2, test.use = "wilcox", min.pct = 0.1, logFC = 0.25) {

	cells.1 <- WhichCells(object = object, ident = ident.1)
	exp1 <- apply(object@data[, cells.1, drop = F], 1, function(x) mean(x = expm1(x = x)))

	cells.2 <- WhichCells(object = object, ident = ident.2)
	exp2 <- apply(object@data[, cells.2, drop = F], 1, function(x) mean(x = expm1(x = x)))

	cluster_markers <- FindMarkers(object = object, ident.1 = ident.1, ident.2 = ident.2, min.pct = min.pct, logfc.threshold = logFC, test.use = test.use)
		
	df <- data.frame(geneID = rownames(cluster_markers),
					 pct.1 = cluster_markers$pct.1, pct.2 = cluster_markers$pct.2,
					 avg_exp.1 = exp1[rownames(cluster_markers)], avg_exp.2 = exp2[rownames(cluster_markers)],
					 avg_log2FC = cluster_markers$avg_logFC/log(2),
					 p_val = cluster_markers$p_val,
					 p_val_adj = cluster_markers$p_val_adj)
	write.table(df, file=out.name, sep="\t", quote=FALSE, row.names = FALSE)

	return()
}


ConvertGeneMarkersTableWithoutExpr <- function(object, old.table, new.table, ident.1, ident.2 = NULL, clusters = NULL) {

	if (length(clusters) == 0) {
		clusters <- levels(object@ident)
	}

	df.old <- read.table(old.table, header=TRUE, sep="\t")

	if (is.null(ident.2)) {
		ident.2 = setdiff(clusters, ident.1)
	}

	cells.1 <- WhichCells(object = object, ident = ident.1)
	exp1 <- apply(object@data[, cells.1, drop = F], 1, function(x) mean(x = expm1(x = x)))

	cells.2 <- WhichCells(object = object, ident = ident.2)
	exp2 <- apply(object@data[, cells.2, drop = F], 1, function(x) mean(x = expm1(x = x)))

	df.new <- data.frame(geneID = rownames(df.old),
					 pct.1 = df.old$pct.1, pct.2 = df.old$pct.2,
					 avg_exp.1 = exp1[rownames(df.old)], avg_exp.2 = exp2[rownames(df.old)],
					 avg_log2FC = df.old$avg_logFC/log(2),
					 p_val = df.old$p_val,
					 p_val_adj = df.old$p_val_adj)
	write.table(df.new, file=new.table, sep="\t", quote=FALSE, row.names = FALSE)

	return()
}


ConvertGeneMarkersTableWithExpr <- function(old.table, new.table) {

	df.old <- read.table(old.table, header=TRUE, sep="\t")

	df.new <- data.frame(geneID = rownames(df.old),
					 pct.1 = df.old$pct.1, pct.2 = df.old$pct.2,
					 avg_exp.1 = df.old$avg_exp.1, avg_exp.2 = df.old$avg_exp.2,
					 avg_log2FC = df.old$avg_logFC/log(2),
					 p_val = df.old$p_val,
					 p_val_adj = df.old$p_val_adj)
	write.table(df.new, file=new.table, sep="\t", quote=FALSE, row.names = FALSE)

	return()
}


# NB: the row data are used automatically in FindMarkers() when a UMI counts-based method is selected!!!
# NEW TABLE FORMAT
ClusterGeneMarkersByPairsNEW <- function(object, out.dir, res, prefix = NULL, clusters = NULL, test.use = "wilcox", min.pct = 0.1, logFC = 0.25) {

	object <- SetClustersIdent(object, res, prefix)

	if (length(clusters) == 0) {
#		clusters <- sort(unique(object@ident)) # all clusters
		clusters <- levels(object@ident)
	}

	for (i in 1:(length(clusters)-1)) {	
		for (j in (i+1):length(clusters)) {
			deg.table <- paste(out.dir, "/DEG_", test.use, "_cl", clusters[i], "-", clusters[j], "_res", res, ".tsv", sep="")
			GeneMarkersTableNEW(object, out.name = deg.table, ident.1 = clusters[i], ident.2 = clusters[j], test.use = test.use, min.pct = min.pct, logFC = logFC)			
		}
	}

	return()
}


# NB: the row data are used automatically in FindMarkers() when a UMI counts-based method is selected!!!
# NEW TABLE FORMAT
ClusterGeneMarkersVsAllNEW <- function(object, out.dir, res, prefix = NULL, clusters = NULL, test.use = "wilcox", min.pct = 0.1, logFC = 0.25) {

	object <- SetClustersIdent(object, res, prefix)

	if (length(clusters) == 0) {
#		clusters <- sort(unique(object@ident)) # all clusters
		clusters <- levels(object@ident)
	}

	for (i in 1:length(clusters)) {
		deg.table <- paste(out.dir, "/DEG_", test.use, "_cl", clusters[i], "-all_res", res, ".tsv", sep="")
		GeneMarkersTableNEW(object, out.name = deg.table, ident.1 = clusters[i], ident.2 = setdiff(clusters, c(clusters[i])), test.use = test.use, min.pct = min.pct, logFC = logFC)	
	}

	return()
}


# select top genes based on adj-pvalue
FilterClusterGeneMarkers <- function(input.table, output.table, logFC.filt = 1, adjpval.filt = 0.1, min.pct1 = 0.1, min.pct2 = 0.1, num = 100) {

	m <- read.table(input.table, sep="\t", header=TRUE)

	# filter out genes based on abs(logFC), adj. pval, and percentage of cells where the gene is detected
	# IMPORTANT!!! Transform to base 2 log
	epsilon <- 0.00000001
	v <- m[which(abs(m$avg_logFC) > logFC.filt*log(2) - epsilon & m$p_val_adj < adjpval.filt & (m$pct.1 > min.pct1 | m$pct.2 > min.pct2)),]

	# exclude undefined p-values and HIV genes
	v <- v[!is.na(v[,1]),]
	markers <- rownames(v)
	v <- v[!(markers %in% HIV_genes),]

	# print the top num DEG, according to the value of ad. p-value
	num <- min(num, nrow(v))
	write.table(v[1:num,], file=output.table,  sep="\t", quote=FALSE)
	
	return()
}


FilterClusterGeneMarkersPairs <- function(object, out.dir, res, prefix = NULL, clusters = NULL, test.use = "wilcox", logFC.filt = 1, adjpval.filt = 0.1, min.pct1 = 0.1, min.pct2 = 0.1, num = 100) {

	object <- SetClustersIdent(object, res, prefix)

	if (length(clusters) == 0) {
#		clusters <- unique(sort(object@ident)) # all clusters
		clusters <- levels(object@ident)
	}

	# N.B.: here the clustering is performed *without* considering HIV genes
	# Nevertheless, HIV genes might be found among the significant DEG
	for (i in 1:(length(clusters)-1)) {
		for (j in (i+1):length(clusters)) {
			deg.table <- paste(out.dir, "/", test.use, "/DEG_", test.use, "_cl", clusters[i], "-", clusters[j], "_res", res, ".tsv", sep="")
			filt.table <- paste(out.dir, "/", test.use, "/DEG_", test.use, "_cl", clusters[i], "-", clusters[j], "_res", res, "-filtered.tsv", sep="")
			FilterClusterGeneMarkers(deg.table, filt.table, logFC.filt = logFC.filt, adjpval.filt = adjpval.filt, min.pct1 = min.pct1, min.pct2 = min.pct2, num = num)
		}
	}
	
	return()
}


FilterClusterGeneMarkersAll <- function(object, out.dir, res, prefix = NULL, clusters = NULL, test.use = "wilcox", logFC.filt = 1, adjpval.filt = 0.1, min.pct1 = 0.1, min.pct2 = 0.1, num = 100) {

	object <- SetClustersIdent(object, res, prefix)

	if (length(clusters) == 0) {
#		clusters <- unique(sort(object@ident)) # all clusters
		clusters <- levels(object@ident)
	}

	# N.B.: here the clustering is performed *without* considering HIV genes
	# Nevertheless, HIV genes might be found among the significant DEG
	for (i in 1:length(clusters)) {
		deg.table <- paste(out.dir, "/", test.use, "/DEG_", test.use, "_cl", clusters[i], "-all_res", res, ".tsv", sep="")
		filt.table <- paste(out.dir, "/", test.use, "/DEG_", test.use, "_cl", clusters[i], "-all_res", res, "-filtered.tsv", sep="")
		FilterClusterGeneMarkers(deg.table, filt.table, logFC.filt = logFC.filt, adjpval.filt = adjpval.filt, min.pct1 = min.pct1, min.pct2 = min.pct2, num = num)
	}
	
	return()
}


FilterDEGtable <- function(table, logFC.filt = 1, abs.logFC = TRUE, adjpval.filt = 0.1, min.pct1 = 0.1, min.pct2 = 0.1, num = 100, sort = FALSE) {

	# filter out genes based on abs(logFC), adj. pval, and percentage of cells where the gene is detected
	# IMPORTANT!!! Transform to base 2 log
	epsilon <- 0.00000001
	if (abs.logFC) {
		v <- table[which(abs(table$avg_log2FC) > logFC.filt - epsilon & table$p_val_adj < adjpval.filt & (table$pct.1 > min.pct1 | table$pct.2 > min.pct2)),]
	} else {
		v <- table[which(table$avg_log2FC > logFC.filt - epsilon & table$p_val_adj < adjpval.filt & (table$pct.1 > min.pct1 | table$pct.2 > min.pct2)),]
	}

	# exclude undefined p-values and HIV genes
	v <- v[!is.na(v[,1]),]
	v <- v[!(v$geneID %in% HIV_genes),]

	# print the top num DEG, according to the value of ad. p-value
	num <- min(num, nrow(v))
	v <- v[1:num,]
	if (sort) {
		v <- v[order(v$avg_log2FC, decreasing=TRUE),]
	}
	
	return(v)
}


SelectUpDown <- function (table) {

	g <- rep("upregulated", nrow(table))
	g[which(table$avg_log2FC < 0)] <- "downregulated"
	df <- data.frame(geneID = table$geneID,	group = g)
	
	return(df)
}


FilterClusterGeneMarkersNEW <- function(input.table, output.table, logFC.filt = 1, abs.logFC = TRUE, adjpval.filt = 0.1, min.pct1 = 0.1, min.pct2 = 0.1, num = 100, sort = FALSE) {

	m <- read.table(input.table, sep="\t", header=TRUE)
	v <- FilterDEGtable(table = m, logFC.filt = logFC.filt, abs.logFC = abs.logFC, adjpval.filt = adjpval.filt, min.pct1 = min.pct1, min.pct2 = min.pct2, num = num, sort = sort)
	write.table(v, file=output.table,  sep="\t", quote=FALSE, row.names = FALSE)
	
	return()
}


FilterClusterGeneMarkersPairsNEW <- function(object, out.dir, res, prefix = NULL, clusters = NULL, test.use = "wilcox", logFC.filt = 1, adjpval.filt = 0.1, min.pct1 = 0.1, min.pct2 = 0.1, num = 100) {

	object <- SetClustersIdent(object, res, prefix)

	if (length(clusters) == 0) {
#		clusters <- sort(unique(object@ident)) # all clusters
		clusters <- levels(object@ident)
	}

	# N.B.: here the clustering is performed *without* considering HIV genes
	# Nevertheless, HIV genes might be found among the significant DEG
	for (i in 1:(length(clusters)-1)) {
		for (j in (i+1):length(clusters)) {
			deg.table <- paste(out.dir, "/", test.use, "/DEG_", test.use, "_cl", clusters[i], "-", clusters[j], "_res", res, ".tsv", sep="")
			filt.table <- paste(out.dir, "/", test.use, "/DEG_", test.use, "_cl", clusters[i], "-", clusters[j], "_res", res, "-filtered.tsv", sep="")
			FilterClusterGeneMarkersNEW(deg.table, filt.table, logFC.filt = logFC.filt, adjpval.filt = adjpval.filt, min.pct1 = min.pct1, min.pct2 = min.pct2, num = num)
		}
	}
	
	return()
}


FilterClusterGeneMarkersAllNEW <- function(object, out.dir, res, prefix = NULL, clusters = NULL, test.use = "wilcox", logFC.filt = 1, adjpval.filt = 0.1, min.pct1 = 0.1, min.pct2 = 0.1, num = 100) {

	object <- SetClustersIdent(object, res, prefix)

	if (length(clusters) == 0) {
#		clusters <- sort(unique(object@ident)) # all clusters
		clusters <- levels(object@ident)
	}

	# N.B.: here the clustering is performed *without* considering HIV genes
	# Nevertheless, HIV genes might be found among the significant DEG
	for (i in 1:length(clusters)) {
		deg.table <- paste(out.dir, "/", test.use, "/DEG_", test.use, "_cl", clusters[i], "-all_res", res, ".tsv", sep="")
		filt.table <- paste(out.dir, "/", test.use, "/DEG_", test.use, "_cl", clusters[i], "-all_res", res, "-filtered.tsv", sep="")
		FilterClusterGeneMarkersNEW(deg.table, filt.table, logFC.filt = logFC.filt, adjpval.filt = adjpval.filt, min.pct1 = min.pct1, min.pct2 = min.pct2, num = num)
	}
	
	return()
}


GOPreprocessing <- function(deg.table, filt.table, logFC.filt = 1, adjpval.filt = 0.1) {

	m <- read.table(deg.table, sep="\t", header=TRUE)
	t <- FilterDEGtable(table = m, logFC.filt = logFC.filt, abs.logFC = FALSE, adjpval.filt = adjpval.filt, num = 1000000)
	write.table(t$geneID, file=filt.table, quote=FALSE, row.names = FALSE, col.names = FALSE)
	
	return()
}


GOPreprocessingClusterByPairs <- function(object, out.dir, res, prefix = NULL, clusters = NULL, test.use = "wilcox", logFC.filt = 1, adjpval.filt = 0.1) {

	object <- SetClustersIdent(object, res, prefix)

	if (length(clusters) == 0) {
		clusters <- levels(object@ident)
	}

	# N.B.: here the clustering is performed *without* considering HIV genes
	# Nevertheless, HIV genes might be found among the significant DEG
	for (i in 1:(length(clusters)-1)) {
		for (j in (i+1):length(clusters)) {
			deg.table <- paste(out.dir, "/", test.use, "/DEG_", test.use, "_cl", clusters[i], "-", clusters[j], "_res", res, ".tsv", sep="")
			filt.table <- paste(out.dir, "/", test.use, "/GOprep/list_DEG_", test.use, "_cl", clusters[i], "-", clusters[j], "_res", res, "-logFC", logFC.filt, "-adjpval", adjpval.filt, ".txt", sep="")
			GOPreprocessing(deg.table, filt.table, logFC.filt, adjpval.filt)
		}
	}
	
	return()
}


GOPreprocessingClusterVsAll <- function(object, out.dir, res, prefix = NULL, clusters = NULL, test.use = "wilcox", logFC.filt = 1, adjpval.filt = 0.1) {

	object <- SetClustersIdent(object, res, prefix)

	if (length(clusters) == 0) {
		clusters <- levels(object@ident)
	}

	# N.B.: here the clustering is performed *without* considering HIV genes
	# Nevertheless, HIV genes might be found among the significant DEG
	for (i in 1:length(clusters)) {
		deg.table <- paste(out.dir, "/", test.use, "/DEG_", test.use, "_cl", clusters[i], "-all_res", res, ".tsv", sep="")
		filt.table <- paste(out.dir, "/", test.use, "/GOprep/list_DEG_", test.use, "_cl", clusters[i], "-all_res", res, "-logFC", logFC.filt, "-adjpval", adjpval.filt, ".txt", sep="")
		GOPreprocessing(deg.table, filt.table, logFC.filt, adjpval.filt)
	}
	
	return()
}


GOenrichPreprocessing <- function(deg.table, filt.table, logFC.filt = 1, adjpval.filt = 0.1) {

	m <- read.table(deg.table, sep="\t", header=TRUE)
	t <- FilterDEGtable(table = m, logFC.filt = logFC.filt, adjpval.filt = adjpval.filt, num = 1000000)
	df <- SelectUpDown(table = t)
	write.table(df, file=filt.table, quote=FALSE, row.names = FALSE, col.names = TRUE, sep="\t")
	
	return()
}


GOenrichPreprocessingClusterByPairs <- function(object, out.dir, res, prefix = NULL, clusters = NULL, test.use = "wilcox", logFC.filt = 1, adjpval.filt = 0.1) {

	object <- SetClustersIdent(object, res, prefix)

	if (length(clusters) == 0) {
		clusters <- levels(object@ident)
	}

	# N.B.: here the clustering is performed *without* considering HIV genes
	# Nevertheless, HIV genes might be found among the significant DEG
	for (i in 1:(length(clusters)-1)) {
		for (j in (i+1):length(clusters)) {
			deg.table <- paste(out.dir, "/", test.use, "/DEG_", test.use, "_cl", clusters[i], "-", clusters[j], "_res", res, ".tsv", sep="")
			filt.table <- paste(out.dir, "/", test.use, "/GOenrichprep/list_DEG_", test.use, "_cl", clusters[i], "-", clusters[j], "_res", res, "-logFC", logFC.filt, "-adjpval", adjpval.filt, ".txt", sep="")
			GOenrichPreprocessing(deg.table, filt.table, logFC.filt, adjpval.filt)
		}
	}
	
	return()
}


GOenrichPreprocessingClusterVsAll <- function(object, out.dir, res, prefix = NULL, clusters = NULL, test.use = "wilcox", logFC.filt = 1, adjpval.filt = 0.1) {

	object <- SetClustersIdent(object, res, prefix)

	if (length(clusters) == 0) {
		clusters <- levels(object@ident)
	}

	# N.B.: here the clustering is performed *without* considering HIV genes
	# Nevertheless, HIV genes might be found among the significant DEG
	for (i in 1:length(clusters)) {
		deg.table <- paste(out.dir, "/", test.use, "/DEG_", test.use, "_cl", clusters[i], "-all_res", res, ".tsv", sep="")
		filt.table <- paste(out.dir, "/", test.use, "/GOenrichprep/list_DEG_", test.use, "_cl", clusters[i], "-all_res", res, "-logFC", logFC.filt, "-adjpval", adjpval.filt, ".txt", sep="")
		GOenrichPreprocessing(deg.table, filt.table, logFC.filt, adjpval.filt)
	}
	
	return()
}


# extract statistics on the DE analysis for a specific list of genes
SelectClusterGeneMarkersAll <- function(object, out.dir, res, prefix = NULL, clusters = NULL, test.use = "wilcox", gene.list, selection.name, plot = FALSE) {

	object <- SetClustersIdent(object, res, prefix)

	if (length(clusters) == 0) {
#		clusters <- unique(sort(object@ident)) # all clusters
		clusters <- levels(object@ident)
	}

	# N.B.: here the clustering is performed *without* considering HIV genes
	# Nevertheless, HIV genes might be found among the significant DEG
	for (i in 1:length(clusters)) {
		deg.table <- paste(out.dir, "/", test.use, "/DEG_", test.use, "_cl", clusters[i], "-all_res", res, ".tsv", sep="")
		filt.table <- paste(out.dir, "/", test.use, "/DEG_", test.use, "_cl", clusters[i], "-all_res", res, "-", selection.name, ".tsv", sep="")
		table <- read.table(deg.table, sep="\t", header = TRUE)
		sel_table <- table[rownames(table) %in% gene.list, ]
		write.table(sel_table, file=filt.table,  sep="\t", quote=FALSE)
		plot.name <- paste(out.dir, "/", test.use, "/heatmap_DEG_", test.use, "_cl", clusters[i], "-all_res", res, "-", selection.name, ".pdf", sep="")
		if (plot) {
			genes <- rownames(sel_table)
			cells <- colnames(object@data)[object@ident %in% sort(clusters)]
			pdf(plot.name, width = 5, height = 4)
			print(DoHeatmap(object = object, genes.use = genes, cells.use = cells, group.order = clusters, rotate.key = TRUE, slim.col.label = TRUE))
			dev.off()
		}
	}
	
	return()
}


HeatmapPlot <- function(object, filt.table, plot.name, cells = NULL) {
		
	table <- read.table(filt.table, sep="\t", header = TRUE)
	genes <- rownames(table)
	if (length(genes) == 0) next

	pdf(plot.name, width = 5, height = 7)	
	print(DoHeatmap(object = object, cells.use = cells, genes.use = genes, rotate.key = TRUE, slim.col.label = TRUE))
	dev.off()

	return()
}


# plot genes on a heatmap
HeatmapPlotGeneList <- function(object, plot.name, clusters = NULL, res, prefix = NULL, genes, sort = FALSE) {

	object <- SetClustersIdent(object, res, prefix)

	if (length(clusters) == 0) {
#		clusters <- unique(sort(object@ident)) # all clusters
		clusters <- levels(object@ident)
	}

	cells <- colnames(object@data)[object@ident %in% sort(clusters)]

	# FIXME is this really sorting by clustering??
	if (sort) {
		heatmap_table <- object@scale.data[rownames(object@scale.data) %in% genes, cells]
		d <- dist(heatmap_table)
		c <- hclust(d)
		genes <- genes[c$order]
	}

	height <- 0.5 + 0.12*length(genes)
	pdf(plot.name, width = 10, height = height)	
	print(DoHeatmap(object = object, cells.use = cells, genes.use = genes, rotate.key = TRUE, slim.col.label = TRUE, group.label.rot = FALSE))
	dev.off()

	return()
}


# plot the DEGs of cl-cl comparison
HeatmapPlotPairs <- function(object, out.dir, clusters = NULL, test, res, prefix = NULL) {

	object <- SetClustersIdent(object, res, prefix)

	if (length(clusters) == 0) {
#		clusters <- unique(sort(object@ident)) # all clusters
		clusters <- levels(object@ident)
	}

	for (i in 1:(length(clusters)-1)) {
		for (j in (i+1):length(clusters)) {
			cells <- colnames(object@data)[object@ident %in% sort(c(clusters[i], clusters[j]))]
			deg.table.filt <- paste(out.dir, "/", test, "/DEG_", test, "_cl", clusters[i], "-", clusters[j], "_res", res, "-filtered.tsv", sep="")

			plot.name <- paste(out.dir, "/", test, "/heatmap_DEG_", test, "_cl", clusters[i], "-", clusters[j], "_res", res, ".pdf", sep="")
			HeatmapPlot(object, deg.table.filt, plot.name, cells = cells)
		}
	}
	return()
}


# plot the DEGs of cl-all comparison (for a specific cluster)
HeatmapPlotAll <- function(object, out.dir, clusters = NULL, test, res, prefix = NULL) {

	object <- SetClustersIdent(object, res, prefix)

	if (length(clusters) == 0) {
#		clusters <- unique(sort(object@ident)) # all clusters
		clusters <- levels(object@ident)
	}

	cells <- colnames(object@data)[object@ident %in% sort(clusters)]
	for (i in 1:length(clusters)) {
		cells <- colnames(object@data)[object@ident %in% sort(clusters)]
		deg.table.filt <- paste(out.dir, "/", test, "/DEG_", test, "_cl", clusters[i], "-all_res", res, "-filtered.tsv", sep="")

		plot.name <- paste(out.dir, "/", test, "/heatmap_DEG_", test, "_cl", clusters[i], "-all_res", res, ".pdf", sep="")
		HeatmapPlot(object, deg.table.filt, plot.name, cells = cells)
	}
	return()
}


HeatmapPlotNEW <- function(object, filt.table, plot.name, cells = NULL, revert = FALSE, width = 5, height = 7) {
		
	table <- read.table(filt.table, sep="\t", header = TRUE)
	genes <- table$geneID
	if (length(genes) == 0) next

	pdf(plot.name, width = width, height = height)	
	print(DoHeatmap(object = object, cells.use = cells, genes.use = genes, rotate.key = TRUE, slim.col.label = TRUE))
	dev.off()

	return()
}


# plot the DEGs of cl-cl comparison
HeatmapPlotPairsNEW <- function(object, out.dir, clusters = NULL, test, res, prefix = NULL) {

	object <- SetClustersIdent(object, res, prefix)

	if (length(clusters) == 0) {
		clusters <- levels(object@ident) # all clusters
	}

	dir.create(paste(out.dir, "/", test, "/heatmap", sep=""), showWarnings = FALSE)

	for (i in 1:(length(clusters)-1)) {
		for (j in (i+1):length(clusters)) {
			cells <- colnames(object@data)[object@ident %in% sort(c(clusters[i], clusters[j]))]
			deg.table.filt <- paste(out.dir, "/", test, "/DEG_", test, "_cl", clusters[i], "-", clusters[j], "_res", res, "-filtered.tsv", sep="")

			plot.name <- paste(out.dir, "/", test, "/heatmap/heatmap_DEG_", test, "_cl", clusters[i], "-", clusters[j], "_res", res, ".pdf", sep="")
			HeatmapPlotNEW(object, deg.table.filt, plot.name, cells = cells)
		}
	}
	return()
}


# plot the DEGs of cl-all comparison (for a specific cluster)
HeatmapPlotAllNEW <- function(object, out.dir, clusters = NULL, test, res, prefix = NULL, width = 5, height = 7) {

	object <- SetClustersIdent(object, res, prefix)

	if (length(clusters) == 0) {
		clusters <- levels(object@ident) # all clusters
	}

	dir.create(paste(out.dir, "/", test, "/heatmap", sep=""), showWarnings = FALSE)

	cells <- colnames(object@data)[object@ident %in% sort(clusters)]
	for (i in 1:length(clusters)) {
		cells <- colnames(object@data)[object@ident %in% sort(clusters)]
		deg.table.filt <- paste(out.dir, "/", test, "/DEG_", test, "_cl", clusters[i], "-all_res", res, "-filtered.tsv", sep="")

		plot.name <- paste(out.dir, "/", test, "/heatmap/heatmap_DEG_", test, "_cl", clusters[i], "-all_res", res, ".pdf", sep="")
		HeatmapPlotNEW(object, deg.table.filt, plot.name, cells = cells, width = width, height = height)
	}
	return()
}


# If "names" is not empty, then do not take into account FC and pval cut-offs
VolcanoPlotGeneral <- function(table, plot.name, title = "", revert = FALSE, names = NULL, log2FCthresh = 0.5, pvalthresh = 0.05) {

	# convert to log2!!!!
	table$avg_logFC <- table$avg_logFC/log(2)
	
	# flip the two conditions
	if (revert) {
		table$avg_logFC <- -table$avg_logFC
	}

	maxFC <- max(table$avg_logFC[which(is.finite(table$avg_logFC))])
	minFC <- min(table$avg_logFC[which(is.finite(table$avg_logFC))])
	
	non_detectable <- which(table$p_val_adj<MINIMUM)
	if (length(non_detectable) > 0) {
		table$p_val_adj[non_detectable] <- rep(MINIMUM,length(non_detectable))
	}

	inf_pos <- which(!is.finite(table$avg_logFC) && table$avg_logFC > 0)
	inf_neg <- which(!is.finite(table$avg_logFC) && table$avg_logFC < 0)
	if (length(inf_pos) > 0) {
		table$avg_logFC[inf_pos] <- rep(minFC,length(inf_pos))
	}	
	if (length(inf_neg) > 0) {
		table$avg_logFC[inf_neg] <- rep(minFC,length(inf_neg))
	}

	df <- data.frame(x = table$avg_logFC, y = -log10(table$p_val_adj), z = rownames(table))
	# "up" or "down" regulated, within the cutoff (otherwise "none")
	df$type <- factor(ifelse(abs(table$avg_logFC) > log2FCthresh & table$p_val_adj < pvalthresh, ifelse(table$avg_logFC > 0, "up", "down"), "none"))
	# labelling: input names if specified, otherwise all "up" and "down" genes
	if (!is.null(names)) {
		df$lab <- factor(ifelse(rownames(table) %in% names, "lab", "notlab"))
	} else {
		df$lab <- factor(ifelse(abs(table$avg_logFC) > log2FCthresh & table$p_val_adj < pvalthresh, "lab", "nolab"))
	}
	
	# colour for the genes: all, up, down, labelled
	col_nsig <- "gray"
	col_up <- "firebrick3"
	col_down <- "blue"
	col_lab <- "black"

	pdf(plot.name, width=8, height=8)

	# to set equal limits on the x and -x axes
	maxabsFC <- max(abs(minFC), abs(maxFC))

	g <- ggplot(data = df, aes(x = x, y = y)) + theme_bw() + xlab(expression(log[2](FC))) + ylab(expression(-log[10](FDR))) + xlim(-maxabsFC, maxabsFC)
	g <- g + geom_point(colour = col_nsig, size = 2) 
	g <- g + geom_point(data = df[which(df$type == "up"),], colour = col_up, size = 2)
	g <- g + geom_point(data = df[which(df$type == "down"),], colour = col_down, size = 2)
	g <- g + geom_point(data = df[which(df$lab == "lab"),], colour = col_lab, size = 2)
	g <- g + geom_text_repel(data = df[which(df$lab == "lab"),], aes(label = z))
	g <- g + ggtitle(label=title) + theme(plot.title=element_text(size=30)) + theme(text=element_text(size=15))
	print(g)
	
	dev.off()

	return()
}


# Filt table can just be a list of gene names
VolcanoPlot <- function(all.table, filt.table, plot.name, title = "", revert = FALSE, genes.use = NULL) {

	all <- read.table(all.table, sep="\t")
	if (!is.null(genes.use)) {
		all <- all[rownames(all) %in% genes.use,]
	}
	filt <- read.table(filt.table, sep="\t")

	# convert to log2!!!!
	all$avg_logFC <- all$avg_logFC/log(2)

	if (revert) {
		all$avg_logFC <- -all$avg_logFC
	}

	maxFC <- max(all$avg_logFC[which(is.finite(all$avg_logFC))])
	minFC <- min(all$avg_logFC[which(is.finite(all$avg_logFC))])
	
	non_detectable <- which(all$p_val_adj<MINIMUM)
	if (length(non_detectable) > 0) {
		all$p_val_adj[non_detectable] <- rep(MINIMUM,length(non_detectable))
	}

	inf_pos <- which(!is.finite(all$avg_logFC) && all$avg_logFC > 0)
	inf_neg <- which(!is.finite(all$avg_logFC) && all$avg_logFC < 0)
	if (length(inf_pos) > 0) {
		all$avg_logFC[inf_pos] <- rep(minFC,length(inf_pos))
	}	
	if (length(inf_neg) > 0) {
		all$avg_logFC[inf_neg] <- rep(minFC,length(inf_neg))
	}

	df <- data.frame(x = all$avg_logFC, y = -log10(all$p_val_adj), z = rownames(all))
	df$type <- factor(ifelse(all$avg_logFC > 0, "up", "down"))
	df$sel <- factor(ifelse(rownames(all) %in% rownames(filt), "sel", "notsel"))
	
	col_nsig <- "gray"
	col_up <- "firebrick3"
	col_down <- "blue"

	pdf(plot.name, width=8, height=8)

	maxabsFC <- max(abs(minFC), abs(maxFC))

	g <- ggplot(data = df, aes(x = x, y = y)) + theme_bw() + xlab(expression(log[2](FC))) + ylab(expression(-log[10](FDR))) + xlim(-maxabsFC, maxabsFC)
	g <- g + geom_point(colour = col_nsig, size = 2) 
	g <- g + geom_point(data = df[which(df$type == "up" & df$sel == "sel"),], colour = col_up, size = 2)
	g <- g + geom_point(data = df[which(df$type == "down"& df$sel == "sel"),], colour = col_down, size = 2)
	g <- g + geom_text_repel(data = df[which(df$sel == "sel"),], aes(label = z))
	g <- g + ggtitle(label=title) + theme(plot.title=element_text(size=30)) + theme(text=element_text(size=15))
	print(g)
	
	dev.off()

	return()
}


# Filt table can just be a list of gene names
VolcanoPlotBulk <- function(all.table, filt.table, plot.name, title = "", revert = FALSE, genes.use = NULL) {

	all <- read.table(all.table, sep="\t")
	if (!is.null(genes.use)) {
		all <- all[rownames(all) %in% genes.use,]
	}
	filt <- read.table(filt.table, sep="\t")

	if (revert) {
		all$logFC <- -all$logFC
	}

	maxFC <- max(all$logFC[which(is.finite(all$logFC))])
	minFC <- min(all$logFC[which(is.finite(all$logFC))])
	
	non_detectable <- which(all$fdr<MINIMUM)
	if (length(non_detectable) > 0) {
		all$fdr[non_detectable] <- rep(MINIMUM,length(non_detectable))
	}

	inf_pos <- which(!is.finite(all$logFC) && all$logFC > 0)
	inf_neg <- which(!is.finite(all$logFC) && all$logFC < 0)
	if (length(inf_pos) > 0) {
		all$logFC[inf_pos] <- rep(minFC,length(inf_pos))
	}	
	if (length(inf_neg) > 0) {
		all$logFC[inf_neg] <- rep(minFC,length(inf_neg))
	}

	df <- data.frame(x = all$logFC, y = -log10(all$fdr), z = rownames(all))
	df$type <- factor(ifelse(all$logFC > 0, "up", "down"))
	df$sel <- factor(ifelse(rownames(all) %in% rownames(filt), "sel", "notsel"))
	
	col_nsig <- "gray"
	col_up <- "firebrick3"
	col_down <- "blue"

	pdf(plot.name, width=8, height=8)

	maxabsFC <- max(abs(minFC), abs(maxFC))

	g <- ggplot(data = df, aes(x = x, y = y)) + theme_bw() + xlab(expression(log[2](FC))) + ylab(expression(-log[10](FDR))) + xlim(-maxabsFC, maxabsFC)
	g <- g + geom_point(colour = col_nsig, size = 2) 
	g <- g + geom_point(data = df[which(df$type == "up" & df$sel == "sel"),], colour = col_up, size = 2)
	g <- g + geom_point(data = df[which(df$type == "down"& df$sel == "sel"),], colour = col_down, size = 2)
	g <- g + geom_text_repel(data = df[which(df$sel == "sel"),], aes(label = z))
	g <- g + ggtitle(label=title) + theme(plot.title=element_text(size=30)) + theme(text=element_text(size=15))
	print(g)
	
	dev.off()

	return()
}


VolcanoPlotPairs <- function(out.dir, test, clusters, res, genes.use = NULL) {

	for (i in 1:(length(clusters)-1)) {
		for (j in (i+1):length(clusters)) {
			deg.table <- paste(out.dir, "/", test, "/DEG_", test, "_cl", clusters[i], "-", clusters[j], "_res", res, ".tsv", sep="")
			deg.table.filt <- paste(out.dir, "/", test, "/DEG_", test, "_cl", clusters[i], "-", clusters[j], "_res", res, "-filtered.tsv", sep="")	

			volcano.plot <- paste(out.dir, "/", test, "/volcano_DEG_", test, "_cl", clusters[i], "-", clusters[j], "_res", res, ".pdf", sep="")
			title <- paste("cl", clusters[i], "-", clusters[j], sep="")
			VolcanoPlot (deg.table, deg.table.filt, volcano.plot, title, genes.use = genes.use)

			volcano_rev.plot <- paste(out.dir, "/", test, "/volcano_DEG_", test, "_cl", clusters[j], "-", clusters[i], "_res", res, ".pdf", sep="")
			title_rev <- paste("cl", clusters[j], "-", clusters[i], sep="")
			VolcanoPlot (deg.table, deg.table.filt, volcano_rev.plot, title_rev, revert = TRUE, genes.use = genes.use)
		}
	}

	return()
}


VolcanoPlotAll <- function(out.dir, test, clusters, res, genes.use = NULL) {

	for (i in 1:length(clusters)) {
		deg.table <- paste(out.dir, "/", test, "/DEG_", test, "_cl", clusters[i], "-all_res", res, ".tsv", sep="")
		deg.table.filt <- paste(out.dir, "/", test, "/DEG_", test, "_cl", clusters[i], "-all_res", res, "-filtered.tsv", sep="")	

		volcano.plot <- paste(out.dir, "/", test, "/volcano_DEG_", test, "_cl", clusters[i], "-all_res", res, ".pdf", sep="")
		title <- paste("cl", clusters[i], "-all", sep="")
		VolcanoPlot (deg.table, deg.table.filt, volcano.plot, title, genes.use = genes.use)
	}

	return()
}


# Filt table can just be a list of gene names
# color.genes: array of genes to be plot in color
# label.genes: array of genes to be labelled
# VolcanoPlotNEW <- function(all.table, filt.table, plot.name, title = "", revert = FALSE, genes.use = NULL) {
VolcanoPlotNEW <- function(all.table, color.genes1 = NULL, label.genes1 = NULL, color.genes2 = NULL, label.genes2 = NULL, plot.name, title = "", 
							group1.name = "upregulated", group2.name = "downregulated", revert = FALSE, genes.use = NULL, colors = c("gray", "firebrick3", "blue"), point.size = 1, highlight = NULL) {

	all <- read.table(all.table, sep="\t", header=TRUE)
	if (!is.null(genes.use)) {
		all <- all[all$geneID %in% genes.use,]
	}

	if (revert) {
		all$avg_log2FC <- -all$avg_log2FC
	}

	maxFC <- max(all$avg_log2FC[which(is.finite(all$avg_log2FC))])
	minFC <- min(all$avg_log2FC[which(is.finite(all$avg_log2FC))])
	
	non_detectable <- which(all$p_val_adj < MINIMUM)
	if (length(non_detectable) > 0) {
		all$p_val_adj[non_detectable] <- rep(MINIMUM, length(non_detectable))
	}

	inf_pos <- which(!is.finite(all$avg_log2FC) && all$avg_log2FC > 0)
	inf_neg <- which(!is.finite(all$avg_log2FC) && all$avg_log2FC < 0)
	if (length(inf_pos) > 0) {
		all$avg_log2FC[inf_pos] <- rep(minFC,length(inf_pos))
	}	
	if (length(inf_neg) > 0) {
		all$avg_log2FC[inf_neg] <- rep(minFC,length(inf_neg))
	}

	df <- data.frame(x = all$avg_log2FC, y = -log10(all$p_val_adj), z = all$geneID)
	df$type <- factor(ifelse(all$geneID %in% color.genes1, group1.name, ifelse(all$geneID %in% color.genes2, group2.name, "nogroup"))) # gene list 1 or 2
	df$type <- factor(df$type, levels = c("nogroup", group1.name, group2.name))
	df$col <- factor(ifelse(all$geneID %in% c(color.genes1, color.genes2), "col", "notcol"))
	df$lab <- factor(ifelse(all$geneID %in% c(label.genes1, label.genes2), "lab", "notlab"))
	df$hlab <- factor(ifelse(all$geneID %in% highlight, "hlab", "nothlab"))

	pdf(plot.name, width=8.3, height=7)

	maxabsFC <- max(abs(minFC), abs(maxFC))

	g <- ggplot(data = df, aes(x = x, y = y, color=type)) + theme_bw() + xlab(expression(log[2](FC))) + ylab(expression(-log[10](FDR))) + xlim(-maxabsFC, maxabsFC)
	g <- g + geom_point(size=point.size)
	g <- g + scale_color_manual(values=colors)
#	g <- g + geom_text_repel(data = df[which(df$lab == "lab"),], aes(label = z), color="black")
	g <- g + geom_text_repel(data = df[which(df$lab == "lab" & df$hlab == "nothlab"),], aes(label = z), color="black")
	g <- g + geom_label_repel(data = df[which(df$lab == "lab" & df$hlab == "hlab"),], aes(label = z), color="forestgreen", fontface = "bold")
	g <- g + ggtitle(label=title) + theme(plot.title=element_text(size=30)) + theme(text=element_text(size=15))
	print(g)
	
	dev.off()

	return()
}


VolcanoPlotPairsNEW <- function(object, out.dir, test, clusters = NULL, res, prefix = NULL, genes.use = NULL, highlight = NULL) {

	object <- SetClustersIdent(object, res, prefix)

	if (length(clusters) == 0) {
		clusters <- levels(object@ident)
	}

	dir.create(paste(out.dir, "/", test, "/volcano", sep=""), showWarnings = FALSE)

	for (i in 1:(length(clusters)-1)) {
		for (j in (i+1):length(clusters)) {
			deg.table <- paste(out.dir, "/", test, "/DEG_", test, "_cl", clusters[i], "-", clusters[j], "_res", res, ".tsv", sep="")
			deg.table.filt <- paste(out.dir, "/", test, "/DEG_", test, "_cl", clusters[i], "-", clusters[j], "_res", res, "-filtered.tsv", sep="")	
			filt <- read.table(deg.table.filt, sep="\t", header=TRUE)
			label.genes1 <- as.character(filt$geneID[filt$avg_log2FC > 0]) # otherwise they are seen as factors!!
			color.genes1 <- as.character(filt$geneID[filt$avg_log2FC > 0])
			label.genes2 <- as.character(filt$geneID[filt$avg_log2FC < 0]) # otherwise they are seen as factors!!
			color.genes2 <- as.character(filt$geneID[filt$avg_log2FC < 0])

			volcano.plot <- paste(out.dir, "/", test, "/volcano/volcano_DEG_", test, "_cl", clusters[i], "-", clusters[j], "_res", res, ".pdf", sep="")
			title <- paste("cl", clusters[i], "-", clusters[j], sep="")
#			VolcanoPlotNEW (deg.table, deg.table.filt, volcano.plot, title, genes.use = genes.use)
			VolcanoPlotNEW (deg.table, label.genes1 = label.genes1, color.genes1 = color.genes1, label.genes2 = label.genes2, color.genes2 = color.genes2, 
				plot.name = volcano.plot, title = title, genes.use = genes.use, point.size = 2, highlight = highlight)

			volcano_rev.plot <- paste(out.dir, "/", test, "/volcano/volcano_DEG_", test, "_cl", clusters[j], "-", clusters[i], "_res", res, ".pdf", sep="")
			title_rev <- paste("cl", clusters[j], "-", clusters[i], sep="")
#			VolcanoPlotNEW (deg.table, deg.table.filt, volcano_rev.plot, title_rev, revert = TRUE, genes.use = genes.use)
			VolcanoPlotNEW (deg.table, label.genes1 = label.genes2, color.genes1 = color.genes2, label.genes2 = label.genes1, color.genes2 = color.genes1, 
				plot.name = volcano_rev.plot, title = title_rev, revert = TRUE, genes.use = genes.use, point.size = 2, highlight = highlight)
		}
	}

	return()
}


VolcanoPlotAllNEW <- function(object, out.dir, test, clusters = NULL, res, prefix = NULL, genes.use = NULL, highlight = NULL) {

	object <- SetClustersIdent(object, res, prefix)

	if (length(clusters) == 0) {
		clusters <- levels(object@ident)
	}

	dir.create(paste(out.dir, "/", test, "/volcano", sep=""), showWarnings = FALSE)

	for (i in 1:length(clusters)) {
		deg.table <- paste(out.dir, "/", test, "/DEG_", test, "_cl", clusters[i], "-all_res", res, ".tsv", sep="")
		deg.table.filt <- paste(out.dir, "/", test, "/DEG_", test, "_cl", clusters[i], "-all_res", res, "-filtered.tsv", sep="")	
		filt <- read.table(deg.table.filt, sep="\t", header=TRUE) 
		label.genes1 <- as.character(filt$geneID[filt$avg_log2FC > 0]) # otherwise they are seen as factors!!
		color.genes1 <- as.character(filt$geneID[filt$avg_log2FC > 0])
		label.genes2 <- as.character(filt$geneID[filt$avg_log2FC < 0]) # otherwise they are seen as factors!!
		color.genes2 <- as.character(filt$geneID[filt$avg_log2FC < 0])

		volcano.plot <- paste(out.dir, "/", test, "/volcano/volcano_DEG_", test, "_cl", clusters[i], "-all_res", res, ".pdf", sep="")
		title <- paste("cl", clusters[i], "-all", sep="")
#		VolcanoPlotNEW (deg.table, deg.table.filt, volcano.plot, title, genes.use = genes.use)
		VolcanoPlotNEW (deg.table, label.genes1 = label.genes1, color.genes1 = color.genes1, label.genes2 = label.genes2, color.genes2 = color.genes2, 
			plot.name = volcano.plot, title = title, genes.use = genes.use, point.size = 2, highlight = highlight)
	}

	return()
}


# return the avg_log2FC for all genes in object@data
# avg_logFC is defined as log(expr1+1) - log(expr2+1), where log is the natural logarithm
Log2FC <- function (object, ident.1, ident.2, assay.type = "RNA") {

	cells.1 <- WhichCells(object = object, ident = ident.1)
	cells.2 <- WhichCells(object = object, ident = ident.2)
	
	data.test <- GetAssayData(object = object, assay.type = assay.type, slot = "data")
    genes.use <- rownames(x = data.test)

	exp1 <- apply(data.test[, cells.1, drop = F], 1, function(x) mean(x = expm1(x = x)))
	exp2 <- apply(data.test[, cells.2, drop = F], 1, function(x) mean(x = expm1(x = x)))
	logFC <- log(exp1+1) - log(exp2+1)
	log2FC <- logFC/log(2)
	
	genes.use <- rownames(x = data.test)
	to.return <- data.frame(id = genes.use, avg_log2FC = log2FC)
	return(to.return)
}


# return the value of the t-statistics for all genes in object@data
Tstat <- function (object, ident.1, ident.2, assay.type = "RNA") {

	cells.1 <- WhichCells(object = object, ident = ident.1)
	cells.2 <- WhichCells(object = object, ident = ident.2)
    data.test <- GetAssayData(object = object, assay.type = assay.type, slot = "data")
    genes.use <- rownames(x = data.test)
    t_stat <- unlist(x = pblapply(X = genes.use, FUN = function(x) {
        t.test(x = data.test[x, cells.1], y = data.test[x, cells.2])$statistic
    }))
    to.return <- data.frame(id = genes.use, t_stat)
    return(to.return)
}


GSEAPreprocessing <- function (object, ident1, ident2, measure, table) {

	if (measure == "t") {
		# WARNING: this is quite slow!
		df <- Tstat(object, ident.1 = ident1, ident.2 = ident2)
		df <- df[!is.na(df$t_stat),]
	} else {
		if (measure == "f") {
			df <- Log2FC(object, ident.1 = ident1, ident.2 = ident2)
		} else {
			cat(paste("Measure", measure, "is not implemented [GSEAPreprocessing]\n"))
			return()
		}
	}
	df <- df[!(df %in% HIV_genes)]
	write.table(df, file=table, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
}


GSEAPreprocessingClusterByPairs <- function (object, out.dir, res, prefix = NULL, clusters = "", measure = "t") {

	object <- SetClustersIdent(object, res, prefix)

	if (length(clusters) == 0) {
		clusters <- levels(object@ident)
	}

	# N.B.: here the clustering is performed *without* considering HIV genes
	# Nevertheless, HIV genes might be found among the significant DEG
	for (i in 1:(length(clusters)-1)) {
		ident1 <- clusters[i]
		for (j in (i+1):length(clusters)) {
			ident2 <- clusters[j]
			if (measure == "t") {
				id <- "by-tstat"
			} else {
				id <- "by-FC"
			}
			table <- paste(out.dir, "/list_", id, "_cl", clusters[i], "-", clusters[j], "_res", res, ".txt", sep="")
			GSEAPreprocessing (object, ident1, ident2, measure, table)
		}
	}
	
	return()
}


GSEAPreprocessingClusterVsAll <- function (object, out.dir, res, prefix = NULL, clusters = NULL, measure = "t") {

	object <- SetClustersIdent(object, res, prefix)

	if (length(clusters) == 0) {
		clusters <- levels(object@ident)
	}

	if (measure == "t") {
		id <- "by-tstat"
	} else {
		id <- "by-FC"
	}

	# N.B.: here the clustering is performed *without* considering HIV genes
	# Nevertheless, HIV genes might be found among the significant DEG
	for (i in 1:length(clusters)) {
		ident1 <- clusters[i]
		ident2 <- setdiff(clusters, c(clusters[i]))
		table <- paste(out.dir, "/list_", id, "_cl", clusters[i], "-all_res", res, ".txt", sep="")
		GSEAPreprocessing (object, ident1, ident2, measure, table)
	}
	
	return()
}


# correlation between the p-values assigned by the 2 methods
DEGCorrelation <- function(table.names, tests) {

	n <- length(table.names)
	cor_matrix <- diag(rep(1,n))
	rownames(cor_matrix) <- colnames(cor_matrix) <- tests
	for (h in 1:(n-1)) {
		t1 <- read.table(table.names[h], header=TRUE, sep="\t")
		names1 <- rownames(t1)
		t1 <- t1[order(names1),]
		for (k in (h+1):n) {
			t2 <- read.table(table.names[k], header=TRUE, sep="\t")
			names2 <- rownames(t2)
			t2 <- t2[order(names2),]
			# compute the correlation among the common DEGs
			# a different FC may result from normalized and raw gene expression
			p1 <- t1[names1 %in% names2,]$p_val_adj
			p2 <- t2[names2 %in% names1,]$p_val_adj
			cor_matrix[k,h] <- cor_matrix[h,k] <- cor(p1,p2)	
		}
	}

	return(cor_matrix)
}


DEGCorrelationPairs <- function(out.dir, tests, clusters, res, tests.name = "", plot = FALSE) {

	if (plot) {
		plot.name <- paste(out.dir, "/cor_DEG_", tests.name, "_pairs_res", res, ".pdf", sep="")
		pdf(plot.name, width=21, height=30)
#		par(mfrow=c(4,3), mar = c(10,10,10,10), xpd=NA)
		par(mfrow=c(4,3), mar = c(15,15,5,5), xpd=NA)
		colors <- colorRampPalette(c("blue", "white", "red"))(n=100)
	}
	for (i in 1:(length(clusters)-1)) {
		for (j in (i+1):length(clusters)) {
			deg.table.names <- paste(out.dir, "/", tests, "/DEG_", tests, "_cl", clusters[i], "-", clusters[j], "_res", res, ".tsv", sep="")
			cor_matrix <- DEGCorrelation (deg.table.names, tests)
			cor.name <- paste(out.dir, "/cor_DEG_", tests.name, "_cl", clusters[i], "-", clusters[j], "_res", res, ".tsv", sep="")
			write.table(cor_matrix, file = cor.name, quote = FALSE, sep="\t")
			if (plot) {
				col <- colors[round(100*(1+min(cor_matrix))/2):100]
				image(as.matrix(cor_matrix), col=col, xaxt="n", yaxt="n")
				axis(1, at=seq(0,1,length.out=ncol(cor_matrix) ), labels= colnames(cor_matrix), las= 2, cex.axis=3)
				axis(2, at=seq(0,1,length.out=ncol(cor_matrix) ), labels= rownames(cor_matrix), las= 2, cex.axis=3)
				title(main=paste("cl", clusters[i], "-", clusters[j], sep=""), cex.main=4)
			}
		}
	}
	if (plot) {
		image(matrix(1:100, nrow=1, ncol=100), col=colors, xaxt="n", yaxt="n")
		axis(2, at=seq(0,1,length.out=2), labels=c(0,1), las= 2, cex.axis=3)
		title(main="Pearson's correlation", cex.main=4)
		dev.off()	
	}
	
	return()
}



DEGCorrelationAll <- function(out.dir, tests, clusters, res, tests.name = "", plot = FALSE) {

	if (plot) {
		plot.name <- paste(out.dir, "/cor_DEG_", tests.name, "_all_res", res, ".pdf", sep="")
		pdf(plot.name, width=21, height=30)
		par(mfrow=c(4,3), mar = c(15,15,5,5), xpd=NA)
		colors <- colorRampPalette(c("blue", "white", "red"))(n=100)
	}
	for (i in 1:length(clusters)) {
		deg.table.names <- paste(out.dir, "/", tests, "/DEG_", tests, "_cl", clusters[i], "-all_res", res, ".tsv", sep="")
		cor_matrix <- DEGCorrelation (deg.table.names, tests)
		cor.name <- paste(out.dir, "/cor_DEG_", tests.name, "_cl", clusters[i], "-all_res", res, ".tsv", sep="")
		write.table(cor_matrix, file = cor.name, quote = FALSE, sep="\t")
		if (plot) {
			col <- colors[round(100*(1+min(cor_matrix))/2):100]
			image(as.matrix(cor_matrix), col=col, xaxt="n", yaxt="n")
			axis(1, at=seq(0,1,length.out=ncol(cor_matrix) ), labels= colnames(cor_matrix), las= 2, cex.axis=3)
			axis(2, at=seq(0,1,length.out=ncol(cor_matrix) ), labels= rownames(cor_matrix), las= 2, cex.axis=3)
			title(main=paste("cl", clusters[i], "-all", sep=""), cex.main=4)
		}
	}
	if (plot) {
		image(matrix(1:100, nrow=1, ncol=100), col=colors, xaxt="n", yaxt="n")
		axis(2, at=seq(0,1,length.out=2), labels=c(0,1), las= 2, cex.axis=3)
		title(main="Pearson's correlation", cex.main=4)
		dev.off()
	}

	return()
}


# compute the size of the intersection 
IntersectionSize <- function(table.names, tests) {

	n <- length(table.names)
	matrix <- v <- matrix(rep(rep(0,n),n), nrow=n, ncol=n) 
	rownames(matrix) <- colnames(matrix) <- tests
	for (h in 1:(n-1)) {
		t1 <- read.table(table.names[h], header=TRUE, sep="\t")
		n1 <- nrow(t1)
		matrix[h,h] <- n1
		for (k in (h+1):n) {
			t2 <- read.table(table.names[k], header=TRUE, sep="\t")
			n2 <- nrow(t2)
			int <- sum(rownames(t1) %in% rownames(t2))
			matrix[h,k] <- matrix[k,h] <- int
			matrix[k,k] <- n2
		}
	}

	return(matrix)
}


# compute the significance of the size of the intersection of the top DEGs
# all.size = size of the whole gene dataset
IntersectionFisher <- function(table.names, tests, all.size) {

	n <- length(table.names)
	pval_matrix <- v <- matrix(rep(rep(1,n),n), nrow=n, ncol=n) # - diag(rep(1,n))
	rownames(pval_matrix) <- colnames(pval_matrix) <- tests
	for (h in 1:(n-1)) {
		t1 <- read.table(table.names[h], header=TRUE, sep="\t")
		n1 <- nrow(t1)
		for (k in (h+1):n) {
			t2 <- read.table(table.names[k], header=TRUE, sep="\t")
			n2 <- nrow(t2)
			int <- sum(rownames(t1) %in% rownames(t2))
			ext <- all.size - (n1 + n2) + 2*int
			contingency_table <- matrix(c(int, n1-int, n2-int, ext), nrow=2, ncol=2)
			res <- fisher.test(contingency_table)
			pval_matrix[h,k] <- pval_matrix[k,h] <- res[[1]]
		}
	}

	return(pval_matrix)
}


# FIXME: the set of DEGs might not be the same if two different methods are used (applied to raw or to normalized UMI counts)
# Either take the intersection of the whole set of DEGs, or take the max/min size of the two
# TODO: it'd be nice to see the intersection genes for each comparison!!!
# DEGIntersectionPairs <- function(out.dir, tests, clusters, tests.name = "", plot = FALSE) {
#
#	if (plot) {
#		plot.name <- paste(out.dir, "/fisher_DEG_", tests.name, "_pairs_res", res, ".pdf", sep="")
#		pdf(plot.name, width=21, height=30)
#		par(mfrow=c(4,3), mar = c(10,10,10,10), xpd=NA)
#		colors <- colorRampPalette(c("blue", "white"))(n=100)
#	}
#	table.all <- paste(out.dir, "/", tests[1], "/DEG_", tests[1], "_cl", clusters[1], "-", clusters[1], "_res", res, ".tsv", sep="")
#	all.size <- nrow(table.all)
#	for (i in 1:(length(clusters)-1)) {
#		for (j in (i+1):length(clusters)) {
#			deg.table.names <- paste(out.dir, "/", tests, "/DEG_", tests, "_cl", clusters[i], "-", clusters[j], "_res", res, "-filtered.tsv", sep="")
#			fisher_matrix <- IntersectionFisher (deg.table.names, tests, all.size = all.size)
#			fisher.name <- paste(out.dir, "/fisher_DEG_", tests.name, "_cl", clusters[i], "-", clusters[j], "_res", res, "-filtered.tsv", sep="")
#			write.table(fisher_matrix, file = fisher.name, quote = FALSE, sep="\t")
#			if (plot) {
#				col <- colors[round(100*min(fisher_matrix)):100]
#				M <- apply(as.matrix(fisher_matrix), 2, log10)
#				image(M, col=col, xaxt="n", yaxt="n")
#				axis(1, at=seq(0,1,length.out=ncol(fisher_matrix) ), labels= colnames(fisher_matrix), las= 2, cex.axis=3)
#				axis(2, at=seq(0,1,length.out=ncol(fisher_matrix) ), labels= rownames(fisher_matrix), las= 2, cex.axis=3)
#				title(main=paste("cl", clusters[i], "-", clusters[j], sep=""), cex.main=4)
#			}
#		}
#	}
#	if (plot) {
#		image(matrix(1:100, nrow=1, ncol=100), col=colors, xaxt="n", yaxt="n")
#		axis(2, at=seq(0,1,length.out=2), labels=c(0,1), las= 2, cex.axis=3)
#		title(main="p-value", cex.main=4)
#		dev.off()	
#	}
#	
#	return()
# }



# FIXME: the set of DEGs might not be the same if two different methods are used (applied to raw or to normalized UMI counts)
# Either take the intersection of the whole set of DEGs, or take the max/min size of the two
# TODO: it'd be nice to see the intersection genes for each comparison!!!
DEGIntersectionPairs <- function(out.dir, tests, clusters, res, tests.name = "", plot = FALSE) {

	colors <- colorRampPalette(c("blue", "white"))(n=100)
	table.all <- paste(out.dir, "/", tests[1], "/DEG_", tests[1], "_cl", clusters[1], "-", clusters[1], "_res", res, ".tsv", sep="")
	all.size <- nrow(table.all)
	for (i in 1:(length(clusters)-1)) {
		for (j in (i+1):length(clusters)) {
			deg.table.names <- paste(out.dir, "/", tests, "/DEG_", tests, "_cl", clusters[i], "-", clusters[j], "_res", res, "-filtered.tsv", sep="")
			intersection_matrix <- IntersectionSize (deg.table.names, tests)
			intersection.name <- paste(out.dir, "/intersect_DEG_", tests.name, "_cl", clusters[i], "-", clusters[j], "_res", res, "-filtered.tsv", sep="")
			write.table(intersection_matrix, file = intersection.name, quote = FALSE, sep="\t")
			fisher_matrix <- IntersectionFisher (deg.table.names, tests, all.size = all.size)
			fisher.name <- paste(out.dir, "/fisher_DEG_", tests.name, "_cl", clusters[i], "-", clusters[j], "_res", res, "-filtered.tsv", sep="")
			write.table(fisher_matrix, file = fisher.name, quote = FALSE, sep="\t")
			if (plot) {
				col <- colors[round(100*min(fisher_matrix)):100]
				plot.name <- paste(out.dir, "/fisher_DEG_", tests.name, "_cl", clusters[i], "-", clusters[j], "_res", res, "-filtered.pdf", sep="")
				M <- apply(as.matrix(fisher_matrix), 2, log10)
				pdf(plot.name, width=6, height=7)
				par(mar=c(12,12,5,1))
				image(M, col=col, xaxt="n", yaxt="n")
				axis(1, at=seq(0,1,length.out=ncol(fisher_matrix) ), labels= colnames(fisher_matrix), las= 2, cex.axis=3)
				axis(2, at=seq(0,1,length.out=ncol(fisher_matrix) ), labels= rownames(fisher_matrix), las= 2, cex.axis=3)
				title(main=paste("cl", clusters[i], "-", clusters[j], sep=""), cex.main=4)
				dev.off()
			}
		}
	}
	if (plot) {
		plot.name <- paste(out.dir, "/fisher_DEG_legend.pdf", sep="")
		pdf(plot.name, width=3, height=5)
		image(matrix(1:100, nrow=1, ncol=100), col=colors, xaxt="n", yaxt="n")
		axis(2, at=seq(0,1,length.out=2), labels=c(0,1), las= 2, cex.axis=3)
		title(main="p-value", cex.main=4)
		dev.off()
	}
	
	return()
}


# DEGIntersectionAll <- function(out.dir, tests, clusters, tests.name = "", plot = FALSE) {
#
#	if (plot) {
#		plot.name <- paste(out.dir, "/fisher_DEG_", tests.name, "_all_res", res, ".pdf", sep="")
#		pdf(plot.name, width=21, height=30)
#		par(mfrow=c(4,3), mar = c(10,10,10,10), xpd=NA)
#		colors <- colorRampPalette(c("blue", "white"))(n=100)
#	}
#	table.all <- paste(out.dir, "/", tests[1], "/DEG_", tests[1], "_cl", clusters[1], "-", clusters[1], "_res", res, ".tsv", sep="")
#	all.size <- nrow(table.all)
#	for (i in 1:length(clusters)) {
#		deg.table.names <- paste(out.dir, "/", tests, "/DEG_", tests, "_cl", clusters[i], "-all_res", res, "-filtered.tsv", sep="")
#		fisher_matrix <- DEGIntersection (deg.table.names, tests, all.size = all.size)
#		fisher.name <- paste(out.dir, "/fisher_DEG_", tests.name, "_cl", clusters[i], "-all_res", res, "-filtered.tsv", sep="")
#		write.table(fisher_matrix, file = fisher.name, quote = FALSE, sep="\t")
#		if (plot) {
#			col <- colors[round(100*min(fisher_matrix)):100]
#			M <- apply(as.matrix(fisher_matrix), 2, log10)
#			image(M, col=col, xaxt="n", yaxt="n")
#			axis(1, at=seq(0,1,length.out=ncol(fisher_matrix) ), labels= colnames(fisher_matrix), las= 2, cex.axis=3)
#			axis(2, at=seq(0,1,length.out=ncol(fisher_matrix) ), labels= rownames(fisher_matrix), las= 2, cex.axis=3)
#			title(main=paste("cl", clusters[i], "-all", sep=""), cex.main=4)
#		}
#	}
#	if (plot) {
#		image(matrix(1:100, nrow=1, ncol=100), col=colors, xaxt="n", yaxt="n")
#		axis(2, at=seq(0,1,length.out=2), labels=c(0,1), las= 2, cex.axis=3)
#		title(main="p-value", cex.main=4)
#		dev.off()	
#	}
#
#	return()
# }


DEGIntersectionAll <- function(out.dir, tests, clusters, res, tests.name = "", plot = FALSE) {

	colors <- colorRampPalette(c("blue", "white"))(n=100)
	table.all <- paste(out.dir, "/", tests[1], "/DEG_", tests[1], "_cl", clusters[1], "-", clusters[1], "_res", res, ".tsv", sep="")
	all.size <- nrow(table.all)
	for (i in 1:length(clusters)) {
		deg.table.names <- paste(out.dir, "/", tests, "/DEG_", tests, "_cl", clusters[i], "-all_res", res, "-filtered.tsv", sep="")
		intersection_matrix <- IntersectionSize (deg.table.names, tests)
		intersection.name <- paste(out.dir, "/intersect_DEG_", tests.name, "_cl", clusters[i], "-all_res", res, "-filtered.tsv", sep="")
		write.table(intersection_matrix, file = intersection.name, quote = FALSE, sep="\t")
		fisher_matrix <- IntersectionFisher (deg.table.names, tests, all.size = all.size)
		fisher.name <- paste(out.dir, "/fisher_DEG_", tests.name, "_cl", clusters[i], "-all_res", res, "-filtered.tsv", sep="")
		write.table(fisher_matrix, file = fisher.name, quote = FALSE, sep="\t")
		if (plot) {
			col <- colors[round(100*min(fisher_matrix)):100]
			plot.name <- paste(out.dir, "/fisher_DEG_", tests.name, "_cl", clusters[i], "-all_res", res, "-filtered.pdf", sep="")
			M <- apply(as.matrix(fisher_matrix), 2, log10)
			pdf(plot.name, width=6, height=7)	
			par(mar=c(12,12,5,1))
			image(M, col=col, xaxt="n", yaxt="n")
			axis(1, at=seq(0,1,length.out=ncol(fisher_matrix) ), labels= colnames(fisher_matrix), las= 2, cex.axis=3)
			axis(2, at=seq(0,1,length.out=ncol(fisher_matrix) ), labels= rownames(fisher_matrix), las= 2, cex.axis=3)
			title(main=paste("cl", clusters[i], "-all", sep=""), cex.main=4)
			dev.off()
		}
	}

	return()
}


# Read UMI count information from a list of barcodes
# File format:
#   CB	UMIcount
# (with header)
# This can be done in two ways:
# - at the sample level, where the barcodes in df are possibly shared between 
#   different samples, and barcode.prefix is added to have unique barcodes
# - at the multi-sample level, where the barcodes already contain the sample
#   information
# In the first case, if the meta.data field is already present, the function 
# sums up the count to the existing count, hence it is assumed that the same 
# sample will not be loaded twice
LoadUMIcount <- function(object, df, meta.name, barcode.prefix = "") {

	if (nrow(df) == 0) {
		return(object)
	}

	# remove the leading -1 at the end (if it exists)
	df$CB <- gsub("-1", "", df$CB)
	if (barcode.prefix != "") {
		df$CB <- paste(barcode.prefix, df$CB, sep="_")
	}

	# init umicounts for the new metadata
	barcodes <- colnames(object@data)
	umicounts <- rep(0, length(barcodes))

	# fill metadata with the input umicounts
	a <- match(df$CB, barcodes)
	if (sum(!is.na(a)) > 0) {
		umicounts[a[!is.na(a)]] <- df$UMIcount[df$CB %in% barcodes]
		if (!(meta.name %in% colnames(object@meta.data))) {
			df_seurat <- data.frame(umicounts)
			rownames(df_seurat) <- barcodes
			colnames(df_seurat) <- meta.name
			object <- AddMetaData(object = object, metadata = df_seurat)
		} else {
			object@meta.data[meta.name] <- object@meta.data[meta.name] + umicounts;
		}
	}

	return(object)
}









