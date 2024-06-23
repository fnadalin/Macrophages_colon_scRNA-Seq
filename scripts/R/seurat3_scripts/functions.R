# Created on 2019/10/11

library("Seurat")
library("future")
library("ggplot2")
library("ggrepel")
library("AUCell")

MINIMUM <- 1e-300



LoadObject <- function(object_file)
{
	l <- load(object_file)
	object <- eval(parse(text=l))
	return(object)
}


# possibly add to an existing object
LoadMultiMatrix <- function(files.list, names.list, out.dir = "", min.cells = 3, min.genes = 200, object = NULL, doublets.list = "") {

	files_list <- as.character(as.array(read.table(files.list)[,1]))
	names_list <- as.character(as.array(read.table(names.list)[,1]))
	doublets_list <- NULL
	if (doublets.list != "") 
		doublets_list <- as.character(as.array(read.table(doublets.list)[,1]))

	if (length(files_list) != length(names_list)) {
		write("<files_list> and <names_list> must have the same length", stderr()) 
		return(1)
	}

	if (!is.null(doublets_list)) {
		if (length(files_list) != length(doublets_list)) {
			write("<files_list> and <doublets_list> must have the same length", stderr()) 
			return(1)
		}
	}

	n <- length(files_list)
	ncells <- 0

	cat("Load the datasets...\n")
	
	# load the dataset
	data <- Read10X(data.dir = files_list[1])
	cat(paste("Matrix dimension:", dim(data), "\n"))
	if (length(doublets_list) > 0) {
		d <- c(as.matrix(read.table(doublets_list[1])))
		data <- data[, which(d == 0)]
		cat(paste("Matrix dimension after doublets removal:", dim(data), "\n"))
	}

	if (n == 1 & is.null(object)) {

		cat(paste("Number of genes: ", nrow(data), "\nNumber of cells: ", ncol(data), "\n", sep=""))
		object <- CreateSeuratObject(counts = data, project = names_list[1], min.cells = min.cells, min.features = min.genes)
		cat(paste("Number of genes after filtering: ", length(rownames(x = object)), "\nNumber of cells after filtering: ", length(colnames(x = object)), "\n", sep=""))

	} else {

		# Create the seurat object
		# Returns a Seurat object with the raw data stored in
		# object@raw.data. object@data, object@meta.data, object@ident, also
		# initialized.
		start <- 1
		if (is.null(object)) {
			object <- CreateSeuratObject(counts = data, project = names_list[1])
			start <- start + 1
		}

		ncells <- ncells + length(colnames(x = object))
		if (n > start) {
			for (i in start:(n-1)) {
				# load the dataset
				data <- Read10X(data = files_list[i])
				cat(paste("Matrix dimension:", dim(data), "\n"))
				if (length(doublets_list) > 0) {
					d <- c(as.matrix(read.table(doublets_list[i])))
					data <- data[, which(d == 0)]
					cat(paste("Matrix dimension after doublets removal:", dim(data), "\n"))
				}
				# Create the seurat object
				obj2 <- CreateSeuratObject(counts = data, project = names_list[i])
				ncells <- ncells + length(colnames(x = object))
				# merge data and filter genes
				# Keep all genes expressed in >= 3 cells (~0.1% of the data). Keep all cells with at least 200 detected genes
				# By doing this at the end, gene filtering is performed based on their expression on the cells of the combined datasets
				if (i == start) {
					object <- merge(x = object, y = obj2, add.cell.ids = c(object@ident, obj2@ident))
				} else {
					object <- merge(x = obj1, y = obj2, add.cell.ids = c(obj1@ident, obj2@ident))
				}
			}
		}
		# load the dataset
		data <- Read10X(data.dir = files_list[n])
		cat(paste("Matrix dimension:", dim(data), "\n"))
		if (length(doublets_list) > 0) {
			d <- c(as.matrix(read.table(doublets_list[n])))
			data <- data[, which(d == 0)]
			cat(paste("Matrix dimension after doublets removal:", dim(data), "\n"))
		}
		# Create the seurat object
		obj2 <- CreateSeuratObject(counts = data, project = names_list[n])
		ncells <- ncells + length(colnames(x = object))

		cat("done\n")

		cat(paste("Number of genes: ", length(rownames(x = object)), "\nNumber of cells: ", ncells, "\n", sep=""))
		cat("Merge into a single object and filter (min.cells = ", min.cells, ", min.genes = ", min.genes, ")\n", sep="")
		object <- merge(x = object, y = obj2, add.cell.ids = c(object@ident, obj2@ident), min.cells = min.cells, min.genes = min.genes)
		cat(paste("Number of genes after filtering: ", length(rownames(x = object)), "\nNumber of cells after filtering: ", length(colnames(x = object)), "\n", sep=""))
	}

	return(object)
}


PlotMetaData <- function (object, violin_plot_features, OUT_DIR, filtered = FALSE, point_size = 1) {

	for (feature in violin_plot_features) {	
		if (!filtered) {
			plot.name <- paste(OUT_DIR, "/", feature, ".pdf", sep="")
		} else {
			plot.name <- paste(OUT_DIR, "/", feature, "_filtered.pdf", sep="")
		}
		pdf(plot.name, width=6, height=6)
		g <- VlnPlot(object = object, features = feature, pt.size = point_size)
		print(g)
		dev.off()
	}

	return()
}


# compute the fraction of transcripts of the input gene group and add this info to a metadata slot
# group_name identifies the gene group in the metadata slot
GeneGroupExpression <- function(object, genes, group_name) {
	
	gene_names <- rownames(x = object)

	M <- GetAssayData(object, slot="counts")

	# FIXME: this does not work if there is just one gene is expressed in some cells!!!
	sum_genes_UMI <- apply(M[gene_names %in% genes,], 2, sum)
#	sum_UMI <- apply(object@raw.data, 2, sum) # N.B. issue with large datasets here. FIXME: use sparse matrix implementation
#	perc_genes_UMI <- sum_genes_UMI/sum_UMI
	perc_genes_UMI <- sum_genes_UMI/object$nCount_RNA # use the UMI computed on ALL the genes!!!!
	log_perc_genes_UMI <- log10(10000*perc_genes_UMI+1)
	
	sum_nGene <- apply(M[gene_names %in% genes,], 2, function(x) sum(x > 0))

	object <- AddMetaData(object = object, metadata = sum_genes_UMI, col.name = paste(group_name, "UMI", sep="."))
	object <- AddMetaData(object = object, metadata = perc_genes_UMI, col.name = paste("perc", group_name, "UMI", sep="."))
	object <- AddMetaData(object = object, metadata = log_perc_genes_UMI, col.name = paste("log.perc", group_name, "UMI", sep="."))
	object <- AddMetaData(object = object, metadata = sum_nGene, col.name = paste(group_name, "genes", sep="."))

	return(object)
}


# FIXME: use mixtools to determine the cut-off between the two modes!!!
BimodalnGeneFilter <- function(object, bimod = TRUE, outliers = FALSE, sd_factor = 2, perc.mito = 0, min.genes = 200, max.genes = 10000, second_max_interval_start = 200, second_max_interval_end_subtract = 1) {

	local_min <- min.genes - 1
	outlier <- max.genes + 1

	cat(paste("Number of cells before nGene filtering: ", length(colnames(x = object)), "\n", sep=""))

	if (bimod) {
		# Filter cells based on sample-specific nGene count
		d <- density(object$nFeature_RNA)
		first_max <- optimize(approxfun(d$x,-d$y), interval=c(min.genes, max.genes))$minimum
		second_max <- optimize(approxfun(d$x,-d$y), interval=c(second_max_interval_start, first_max-second_max_interval_end_subtract))$minimum
		# local minimum between the two modes
		local_min <- optimize(approxfun(d$x,d$y), interval=c(min(first_max,second_max), max(first_max,second_max)))$minimum
	}

	if (outliers) {
		mean_nGene <- mean(object$nFeature_RNA)
		sd_nGene <- sd(object$nFeature_RNA)
		outlier <- mean_nGene + sd_factor*sd_nGene
	}

#	object <- subset(x = object, subset = nFeature_RNA > localmin & nFeature_RNA < outlier)
#	cells.use <- WhichCells(object = object, expression = nFeature_RNA > localmin & nFeature_RNA < outlier)
	object <- SubsetData(object = object, subset.name = "nFeature_RNA", low.threshold = local_min, high.threshold = outlier)
	cat(paste("Number of cells after nGene filtering: ", length(colnames(x = object)), "\n", sep=""))

	if (perc.mito > 0) {
		# Further filter cells on fraction of mito genes 
		object <- subset(x = object, subset = perc.mito.UMI < perc.mito)
		cat(paste("Number of cells after filtering by mitochondrial gene expression: ", length(colnames(x = object)), "\n", sep=""))
	}

	return(object)
}


NormScaleMatrix <- function(object, model = "linear", vars.to.regress = "nCount_RNA", scale = TRUE) { 

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
	if (scale) 
		obj <- ScaleData(object = obj, vars.to.regress = vars.to.regress, model.use=model)

	cat("done\n")

	return(obj)
}


# remove cells with high expression of either one of the genes in the input array
FilterCellsByHighGeneExpression <- function(object, genes, sd_factor = 3) {

	M <- GetAssayData(object, slot="data")

	# compute the mean and sdev at the beginning
	high_thresholds <- low_thresholds <- rep(0,length(genes))
	for (i in 1:length(genes)) {
		if (genes[i] %in% rownames(x = object)) {
			m <- mean(M[genes[i],])
			s <- sd(M[genes[i],])
			high_thresholds[i] <- m + sd_factor*s	
		}
	}	

	for (i in 1:length(genes)) {
		if (genes[i] %in% rownames(x = object)) {
			object <- subset(x = object, subset = genes[i] > low_thresholds[i] & genes[i] < high_thresholds[i])
			cat(paste("Number of cells after ", genes[i], " filtering: ", length(colnames(x = object)), "\n", sep=""))
		}
	}
	
	return(object)
}


ComputeAucell <- function(object, genes, aucMaxRankPerc = 0.2, plot = "aucell.pdf") {

	exprMatrix <- as.matrix(GetAssayData(object, slot="data"))
	pdf(plot)
	cells_rankings <- AUCell_buildRankings(exprMatrix)
	cells_AUC <- AUCell_calcAUC(genes, cells_rankings, aucMaxRank=nrow(cells_rankings)*aucMaxRankPerc)
	dev.off()

	M <- as.matrix(getAUC(cells_AUC))
	for (i in 1:nrow(M)) {
		meta <- data.frame(M[i,])	
		colnames(meta) <- rownames(getAUC(cells_AUC))[i]
		rownames(meta) <- colnames(getAUC(cells_AUC))
		object <- AddMetaData(object = object, metadata = meta)
	}

	return(object)
}


FilterCellsByGeneAuc <- function(object, genes, aucMaxRankPerc = 0.2, auc.cutoff = 0.25) {

	object <- ComputeAucell(object, genes, aucMaxRankPerc)
	object <- subset(x = object, subset = names(genes) < rep(auc.cutoff, length(genes)))
	cat(paste("Number of cells after AUC filtering: ", length(colnames(x = object)), "\n", sep=""))

	return(object)
}


# remove the cells that belong to the specified clusters
# clusters are specified as a comma-separated string
# return the list of selected cells 
# select = TRUE will select the specified identities instead of removing them
FilterOutClusters <- function (object, id, cl_filter, select = FALSE) {

	if (cl_filter == "") {
		return(colnames(GetAssayData(object)))
	}

	cl <- unlist(strsplit(cl_filter, split=","))	

	if (select) {
		filtered_cells <- names(object@active.ident[object@meta.data[[id]] %in% cl])
	} else {
		filtered_cells <- names(object@active.ident[!(object@meta.data[[id]] %in% cl)])
	}

	return(filtered_cells)
}


# reassign cluster names
# renaming is specified as a comma-separated string and refers to clusters 0,1,2... in this order
# specify the array of clusters to be masked (i.e. exlcuded from the renaming
RenameClusters <- function (object, renaming, id, cl_filter = NULL) {

	old.cl.names <- orig.cl.names <- as.character(sort(as.numeric(unique(object@meta.data[[id]]))))
	if (!is.null(cl_filter)) {
		cl <- unlist(strsplit(cl_filter, split=","))
		old.cl.names <- orig.cl.names[!(orig.cl.names %in% cl)]
		object@meta.data[[id]] <- plyr::mapvalues(x = object@meta.data[[id]], from = cl, to = rep("X", length(cl))) # assign X to removed clusters
	}
	new.cl.names <- unlist(strsplit(renaming, split=","))
	object@meta.data[[id]] <- plyr::mapvalues(x = object@meta.data[[id]], from = old.cl.names, to = new.cl.names) # rename only filtered clusters according to the mapping

	return(object)
}


# the logFC is defined as log(expr1+1) - log(expr2+1), where log is the natural logarithm
GeneMarkersTable <- function(object, out.name, ident.1, ident.2, test.use = "wilcox", min.pct = 0.1, logFC = 0.25, cons.by.group = FALSE, check_sanity = TRUE) {

	cells.1 <- WhichCells(object = object, idents = ident.1)
	exp1 <- apply(GetAssayData(object = object)[, cells.1, drop = F], 1, function(x) mean(x = expm1(x = x)))

	cells.2 <- WhichCells(object = object, idents = ident.2)
	exp2 <- apply(GetAssayData(object = object)[, cells.2, drop = F], 1, function(x) mean(x = expm1(x = x)))

	if (!cons.by.group) {
		cluster_markers <- FindMarkers(object = object, ident.1 = ident.1, ident.2 = ident.2, min.pct = min.pct, logfc.threshold = logFC, test.use = test.use)
		df <- data.frame(geneID = rownames(cluster_markers),
						 pct.1 = cluster_markers$pct.1, pct.2 = cluster_markers$pct.2,
						 avg_exp.1 = exp1[rownames(cluster_markers)], avg_exp.2 = exp2[rownames(cluster_markers)],
						 avg_log2FC = cluster_markers$avg_logFC/log(2),
						 p_val = cluster_markers$p_val,
						 p_val_adj = cluster_markers$p_val_adj)
	} else {
		cluster_markers <- FindConservedMarkers(object = object, ident.1 = ident.1, ident.2 = ident.2, min.pct = min.pct, logfc.threshold = logFC, test.use = test.use, grouping.var = "orig.ident")
		df <- data.frame(geneID = rownames(cluster_markers), cluster_markers)
	}
		
	write.table(df, file=out.name, sep="\t", quote=FALSE, row.names = FALSE)

	return()
}


# the logFC is defined as log(expr1+1) - log(expr2+1), where log is the natural logarithm
GeneMarkersTableNEW <- function(object, out.name, ident.1, ident.2, test.use = "wilcox", min.pct = 0.1, logFC = 0.25, cons.by.group = FALSE, check_sanity = TRUE) {

	cells.1 <- WhichCells(object = object, idents = ident.1)
	exp1 <- apply(GetAssayData(object = object, assay = "RNA")[, cells.1, drop = F], 1, function(x) mean(x = expm1(x = x)))

	cells.2 <- WhichCells(object = object, idents = ident.2)
	exp2 <- apply(GetAssayData(object = object, assay = "RNA")[, cells.2, drop = F], 1, function(x) mean(x = expm1(x = x)))

	if (!cons.by.group) {
		# NEW: specify the assay explicitely!!!
		# this is to make sure that the integrated data will NOT be used for differential expression
		cluster_markers <- FindMarkers(object = object, assay = "RNA", ident.1 = ident.1, ident.2 = ident.2, min.pct = min.pct, logfc.threshold = logFC, test.use = test.use)

		df <- data.frame(geneID = rownames(cluster_markers),
						 pct.1 = cluster_markers$pct.1, pct.2 = cluster_markers$pct.2,
						 avg_exp.1 = exp1[rownames(cluster_markers)], avg_exp.2 = exp2[rownames(cluster_markers)],
						 avg_log2FC = cluster_markers$avg_logFC/log(2),
						 p_val = cluster_markers$p_val,
						 p_val_adj = cluster_markers$p_val_adj)
	} else {
		cluster_markers <- FindConservedMarkers(object = object, assay = "RNA", ident.1 = ident.1, ident.2 = ident.2, min.pct = min.pct, logfc.threshold = logFC, test.use = test.use, grouping.var = "orig.ident")
		df <- data.frame(geneID = rownames(cluster_markers), cluster_markers)
	}
		
	write.table(df, file=out.name, sep="\t", quote=FALSE, row.names = FALSE)

	return()
}


# NB: the row data are used automatically in FindMarkers() when a UMI counts-based method is selected!!!
# NEW TABLE FORMAT
ClusterGeneMarkersByPairs <- function(object, out.dir, id = "", clusters = NULL, test.use = "wilcox", min.pct = 0.1, logFC = 0.25, cons.by.group = FALSE) {

	if (id != "") Idents(object) <- id
	if (length(clusters) == 0) clusters <- levels(object@active.ident)

	for (i in 1:(length(clusters)-1)) {	
		for (j in (i+1):length(clusters)) {
			deg.table <- paste(out.dir, "/DEG_", test.use, "_cl", clusters[i], "-", clusters[j], ".tsv", sep="")
#			GeneMarkersTable(object, out.name = deg.table, ident.1 = clusters[i], ident.2 = clusters[j], test.use = test.use, min.pct = min.pct, logFC = logFC, cons.by.group = cons.by.group)	
                        GeneMarkersTableNEW(object, out.name = deg.table, ident.1 = clusters[i], ident.2 = clusters[j], test.use = test.use, min.pct = min.pct, logFC = logFC, cons.by.group = cons.by.group)		
		}
	}

	return()
}


# NB: the row data are used automatically in FindMarkers() when a UMI counts-based method is selected!!!
# NEW TABLE FORMAT
ClusterGeneMarkersVsAll <- function(object, out.dir, id = "", clusters = NULL, test.use = "wilcox", min.pct = 0.1, logFC = 0.25, cons.by.group = FALSE) {

	if (id != "") Idents(object) <- id
	if (length(clusters) == 0) clusters <- levels(object@active.ident)

	for (i in 1:length(clusters)) {
		deg.table <- paste(out.dir, "/DEG_", test.use, "_cl", clusters[i], "-all.tsv", sep="")
#		GeneMarkersTable(object, out.name = deg.table, ident.1 = clusters[i], ident.2 = setdiff(clusters, c(clusters[i])), test.use = test.use, min.pct = min.pct, logFC = logFC, cons.by.group = cons.by.group)
                GeneMarkersTableNEW(object, out.name = deg.table, ident.1 = clusters[i], ident.2 = setdiff(clusters, c(clusters[i])), test.use = test.use, min.pct = min.pct, logFC = logFC, cons.by.group = cons.by.group)	
	}

	return()
}


FilterDEGtable <- function(table, logFC.filt = 1, adjpval.filt = 0.1, min.pct = 0.1, num = 100, sort = FALSE, genes.use = NULL, genes.filter = NULL) {

	# filter out genes based on abs(logFC), adj. pval, and percentage of cells where the gene is detected
	# IMPORTANT!!! Transform to base 2 log
	epsilon <- 0.00000001
	v <- table[which(abs(table$avg_log2FC) > logFC.filt - epsilon & table$p_val_adj < adjpval.filt & (table$pct.1 > min.pct | table$pct.2 > min.pct)),]

	# exclude undefined p-values and HIV genes
	v <- v[!is.na(v[,1]),]

	# exclude gene.filter genes / keep only specified genes
	if (!is.null(genes.filter)) v <- v[!(v$geneID %in% genes.filter),]
	if (!is.null(genes.use)) v <- v[v$geneID %in% genes.use,]

	# print the top num DEG, according to the value of ad. p-value
	num <- min(num, nrow(v))
	v <- v[1:num,]
	if (sort) {
		v <- v[order(v$avg_log2FC, decreasing=TRUE),]
	}
	
	return(v)
}


FilterClusterGeneMarkers <- function(input.table, output.table, logFC.filt = 1, adjpval.filt = 0.1, min.pct = 0.1, num = 100, sort = FALSE, genes.use = NULL, genes.filter = NULL) {

	m <- read.table(input.table, sep="\t", header=TRUE)
	v <- FilterDEGtable(table = m, logFC.filt = logFC.filt, adjpval.filt = adjpval.filt, min.pct = min.pct, num = num, sort = sort, genes.use = genes.use, genes.filter = NULL)
	write.table(v, file=output.table,  sep="\t", quote=FALSE, row.names = FALSE)
	
	return()
}


HeatmapPlot <- function(object, filt.table, plot.name, cells = NULL, width = 5, height = 7) {
		
	table <- read.table(filt.table, sep="\t", header = TRUE)
	genes <- as.character(table$geneID)
	if (length(genes) == 0) next

	pdf(plot.name, width = width, height = height)	
	print(DoHeatmap(object = object, cells = cells, features = genes, label = TRUE, size = 2))
	dev.off()

	return()
}


HeatmapPlotNEW <- function(object, filt.table, plot.name, cells = NULL, width = 5, height = 7) {
		
	table <- read.table(filt.table, sep="\t", header = TRUE)
	genes <- as.character(table$geneID)
	if (length(genes) == 0) next

	pdf(plot.name, width = width, height = height)	
	print(DoHeatmap(object = object, cells = cells, features = genes, label = TRUE, size = 2, assay = "RNA"))
	dev.off()

	return()
}


# Filt table can just be a list of gene names
# color.genes: array of genes to be plot in color
# label.genes: array of genes to be labelled
# VolcanoPlotNEW <- function(all.table, filt.table, plot.name, title = "", revert = FALSE, genes.use = NULL) {
VolcanoPlot <- function(all.table, color.genes = NULL, label.genes = NULL, plot.name, title = "", subtitle = "", group.name = "top50sign",
						revert = FALSE, genes.use = NULL, genes.filter = NULL, colors = c("gray", "blue"), point.size = 1) {

	all <- read.table(all.table, sep="\t", header=TRUE)
	if (!is.null(genes.use)) all <- all[all$geneID %in% genes.use,]
	if (!is.null(genes.filter)) all <- all[!(all$geneID %in% genes.filter),]

	if (revert) all$avg_log2FC <- -all$avg_log2FC

	maxFC <- max(all$avg_log2FC[which(is.finite(all$avg_log2FC))])
	minFC <- min(all$avg_log2FC[which(is.finite(all$avg_log2FC))])
	
	non_detectable <- which(all$p_val_adj < MINIMUM)
	if (length(non_detectable) > 0) 
		all$p_val_adj[non_detectable] <- rep(MINIMUM, length(non_detectable))

	inf_pos <- which(!is.finite(all$avg_log2FC) && all$avg_log2FC > 0)
	inf_neg <- which(!is.finite(all$avg_log2FC) && all$avg_log2FC < 0)
	if (length(inf_pos) > 0) 
		all$avg_log2FC[inf_pos] <- rep(minFC,length(inf_pos))
	if (length(inf_neg) > 0) 
		all$avg_log2FC[inf_neg] <- rep(minFC,length(inf_neg))

	df <- data.frame(x = all$avg_log2FC, y = -log10(all$p_val_adj), z = all$geneID)
	df$type <- factor(ifelse(all$geneID %in% color.genes, group.name, "nogroup"))
	df$type <- factor(df$type, levels = c("nogroup", group.name))
	df$col <- factor(ifelse(all$geneID %in% color.genes, "col", "notcol"))
	df$lab <- factor(ifelse(all$geneID %in% label.genes, "lab", "notlab"))

	pdf(plot.name, width=6.5, height=7)

	maxabsFC <- max(abs(minFC), abs(maxFC))

	g <- ggplot(data = df, aes(x = x, y = y, color=type)) + theme_bw() + xlab(expression(log[2](FC))) + ylab(expression(-log[10](FDR))) + xlim(-maxabsFC, maxabsFC)
	g <- g + geom_point(size=point.size)
	g <- g + scale_color_manual(values=colors)
	g <- g + geom_text_repel(data = df[which(df$lab == "lab"),], aes(label = z), color="black")
	g <- g + ggtitle(label=title, subtitle = subtitle) + theme(plot.title=element_text(size=20), plot.subtitle=element_text(size=10))
	g <- g + theme(text=element_text(size=15)) + theme(legend.position = "none")
	print(g)
	
	dev.off()

	return()
}


VolcanoPlotFilter <- function(all.table, filt.table, plot.name, title = "", subtitle = "", genes.use = NULL, genes.filter = NULL, revert = FALSE) {

	filt <- read.table(filt.table, sep="\t", header=TRUE)

	label.genes <- color.genes <- as.character(filt$geneID) # otherwise they are seen as factors!!

	VolcanoPlot(all.table = all.table, color.genes = color.genes, label.genes = label.genes, 
		plot.name = plot.name, title = title, subtitle = subtitle, genes.use = genes.use, genes.filter = genes.filter, revert = revert)

	return()
}


FilterClusterGeneMarkersPairs <- function(object, out.dir, id = "", clusters = NULL, test.use = "wilcox", logFC.filt = 1, adjpval.filt = 0.1, min.pct = 0.1, num = 100, genes.use = NULL, heatmap = FALSE, volcano = FALSE) {

	if (id != "") Idents(object) <- id
	if (length(clusters) == 0) clusters <- levels(Idents(object))

	dir.create(paste(out.dir, "heatmap", sep="/"), showWarnings = FALSE)
	dir.create(paste(out.dir, "volcano", sep="/"), showWarnings = FALSE)

	for (i in 1:(length(clusters)-1)) {
		for (j in (i+1):length(clusters)) {
			deg.table <- paste(out.dir, "/DEG_", test.use, "_cl", clusters[i], "-", clusters[j], ".tsv", sep="")
			filt.table <- paste(out.dir, "/DEG_", test.use, "_cl", clusters[i], "-", clusters[j], "-filtered.tsv", sep="")
			if (file.exists(deg.table)) {
				FilterClusterGeneMarkers(deg.table, filt.table, logFC.filt = logFC.filt, adjpval.filt = adjpval.filt, min.pct = min.pct, num = num, genes.use = genes.use)
				if (heatmap) {
					cells <- colnames(GetAssayData(object))[Idents(object) %in% sort(c(clusters[i], clusters[j]))]
					plot.name <- paste(out.dir, "/heatmap/heatmap_DEG_", test.use, "_cl", clusters[i], "-", clusters[j], ".pdf", sep="")
					HeatmapPlot(object, filt.table, plot.name, cells = cells)
				}
				if (volcano) {
					volcano.plot <- paste(out.dir, "/volcano/volcano_DEG_", test.use, "_cl", clusters[i], "-", clusters[j], ".pdf", sep="")
					title <- paste("cl", clusters[i], "-", clusters[j], sep="")
					VolcanoPlotFilter(deg.table, filt.table, plot.name = volcano.plot, title = title, genes.use = genes.use)

					volcano_rev.plot <- paste(out.dir, "/volcano/volcano_DEG_", test.use, "_cl", clusters[j], "-", clusters[i], ".pdf", sep="")
					title_rev <- paste("cl", clusters[j], "-", clusters[i], sep="")
					VolcanoPlotFilter(deg.table, filt.table, plot.name = volcano_rev.plot, title = title_rev, revert = TRUE, genes.use = genes.use)
				}
			}
		}
	}
	
	return()
}

FilterClusterGeneMarkersAll <- function(object, out.dir, id = "", clusters = NULL, test.use = "wilcox", logFC.filt = 1, adjpval.filt = 0.1, min.pct = 0.1, num = 100, genes.use = NULL, heatmap = FALSE, volcano = FALSE) {

	dir.create(paste(out.dir, "heatmap", sep="/"))
	dir.create(paste(out.dir, "volcano", sep="/"))

	if (id != "") Idents(object) <- id
	if (length(clusters) == 0) clusters <- levels(Idents(object))
	cells <- colnames(GetAssayData(object))[Idents(object) %in% clusters]

	for (i in 1:length(clusters)) {
		deg.table <- paste(out.dir, "/DEG_", test.use, "_cl", clusters[i], "-all.tsv", sep="")
		filt.table <- paste(out.dir, "/DEG_", test.use, "_cl", clusters[i], "-all-filtered.tsv", sep="")
		if (file.exists(deg.table)) {
		FilterClusterGeneMarkers(deg.table, filt.table, logFC.filt = logFC.filt, adjpval.filt = adjpval.filt, min.pct = min.pct, num = num, genes.use = genes.use)
			if (heatmap) {
				plot.name <- paste(out.dir, "/heatmap/heatmap_DEG_", test.use, "_cl", clusters[i], "-all.pdf", sep="")
				HeatmapPlot(object, filt.table, plot.name, cells = cells, width = 9, height = 7)
			}
			if (volcano) {
				volcano.plot <- paste(out.dir, "/volcano/volcano_DEG_", test.use, "_cl", clusters[i], "-all.pdf", sep="")
				title <- paste("cl", clusters[i], "-all", sep="")
				VolcanoPlotFilter(deg.table, filt.table, plot.name = volcano.plot, title = title, genes.use = genes.use)
			}
		}
	}
	
	return()
}


SamplesComposition <- function(object, out.dir, id_slot = "orig.ident", prefix = NULL, sample.names = NULL, cl.names = NULL) {

	orig_ident <- Idents(object)
	object <- SetIdent(object, value = id_slot)

	if (is.null(sample.names)) 
		sample.names <- unique(orig_ident)

	if (is.null(cl.names)) 
		cl.names <- unique(Idents(object))

	cl_distribution <- matrix(NA, nrow = length(sample.names)*length(cl.names), ncol = 3)
	cl_distribution[,3] <- 0 

	n <- 1
	for (s in sample.names[length(sample.names)-0:(length(sample.names)-1)]) {
		for (c in cl.names) {
			count <- sum(orig_ident == s & Idents(object) == c)
			cl_distribution[n,] <- c(s, c, count)
			n <- n + 1
		}
	}

	cl_distribution <- as.data.frame(cl_distribution)
	colnames(cl_distribution) <- c("sample", "cluster", "num") 
	write.table(cl_distribution, paste(out.dir, "/samples_clusters_composition.tsv", sep=""), sep="\t", row.names=FALSE, quote=FALSE)

	return()
}


BarplotSampleComposition <- function(out.dir, ordering = NULL) {

	# plot cluster distribution
	data <- read.table(paste(out.dir, "/samples_clusters_composition.tsv", sep=""), sep="\t", header=TRUE)
	data$cluster <- as.factor(data$cluster)
	
	if (length(ordering) > 0) {
		data$sample <- factor(data$sample, levels = ordering)
	}

	pdf(paste(out.dir, "/samples_clusters_composition.pdf", sep=""), width=6, height=3)
	print(ggplot(data, aes(x=sample, y=num, fill=cluster)) + geom_bar(stat="identity") + theme_minimal() + coord_flip() + xlab("") + ylab("Number of cells") + theme(axis.text.y=element_text(size=15)))
	dev.off()

	pdf(paste(out.dir, "/samples_clusters_composition_norm.pdf", sep=""), width=6, height=3)
	print(ggplot(data, aes(x=sample, y=num, fill=cluster)) + geom_bar(stat="identity", position="fill") + theme_minimal() + coord_flip() + xlab("") + ylab("Fraction of cells") + theme(axis.text.y=element_text(size=15)))
	dev.off()
	
	return()
}


BarplotClusterComposition <- function(out.dir, ordering = NULL) {

	# plot cluster distribution
	data <- read.table(paste(out.dir, "/samples_clusters_composition.tsv", sep=""), sep="\t", header=TRUE)
	data$sample <- as.factor(data$sample)
	data$cluster <- as.factor(data$cluster)
	
	if (length(ordering) > 0) {
		data$cluster <- factor(data$cluster, levels = ordering)
	}

	pdf(paste(out.dir, "/clusters_samples_composition.pdf", sep=""), 6, height=0.4*length(unique(data$cluster)))
	print(ggplot(data, aes(x=cluster, y=num, fill=sample)) + geom_bar(stat="identity") + theme_minimal() + coord_flip() + xlab("") + ylab("Number of cells") + theme(axis.text.y=element_text(size=15)))
	dev.off()

	pdf(paste(out.dir, "/clusters_samples_composition_norm.pdf", sep=""), width=6, height=0.4*length(unique(data$cluster)))
	print(ggplot(data, aes(x=cluster, y=num, fill=sample)) + geom_bar(stat="identity", position="fill") + theme_minimal() + coord_flip() + xlab("") + ylab("Fraction of cells") + theme(axis.text.y=element_text(size=15)))
	dev.off()
	
	return()
}



