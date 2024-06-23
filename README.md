# Macrophages_colon_scRNA-Seq

This repository contains the code to reproduce the analysis reported in the article:

Chikina *et al.* Macrophages Maintain Epithelium Integrity by Limiting Fungal Product Absorption. *Cell*, 138(2):411-428.e16, 2020. https://doi.org/10.1016/j.cell.2020.08.048.

[System requirements](#system-requirements)  
[Installation](#installation)  
[Instructions for reproducing the analysis](#instructions-for-reproducing-the-analysis)  
[Contact](#contact) 

## System requirements

Upstream cell calling was performed with CellRanger v2.1.1.
The analysis is primarily based on Seurat (v >= 3.0.0).
The code was executed on a GNU/Linux machine with 64GB RAM.

## Instructions for reproducing the analysis

### 1. Data download

The feature-barcode matrices should be downloaded from [GSE146131](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146131) (RAW files) and stored as shown below:

```
$ cd exp/Macrophages_DC_colon_scRNASeq/data/
$ ls Macro_distal
barcodes.tsv.gz  features.tsv.gz  matrix.mtx.gz 
```

The name of the folders in exp/Macrophages_DC_colon_scRNASeq/data/ should correspond to the names reported in exp/Macrophages_DC_colon_scRNASeq/analysis/Seurat/files_list_Macro.
Also, download the mitochondrial and ribosomal gene names from [mito.carta](https://www.broadinstitute.org/files/shared/metabolism/mitocarta/mouse.mitocarta.2.0.html) and [miyazaki](http://ribosome.med.miyazaki-u.ac.jp/):

```
$ cd data/
$ ls
mouse_mitochondrial_genes.txt  mouse_ribosomal_genes.txt
```

### 2. Analysis

The analysis consists of 4 steps:

```
$ cd exp/Macrophages_DC_colon_scRNASeq/analysis/Seurat
$ bash COMMAND_Macro_separate_replicates_SEURAT3
$ bash COMMAND_Macro_merged_replicates_test5_SEURAT3
$ bash COMMAND_optimal_clustering_Macro_SEURAT3
$ bash COMMAND_optimal_clustering_Macro_DEG_SEURAT3
```

each of which is briefly described in the accompanying README file.

## Contact

For feedback or questions about this repository, please contact [Francesca Nadalin](mailto:francesca@ebi.ac.uk). 
