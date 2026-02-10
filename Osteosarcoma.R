##################################################
# STEP 0: Load required libraries
# We need Seurat for analysis and tidyverse ( or dplyr, ggplot2) for data handling & plots
##################################################

library(Seurat)
library(tidyverse)

##################################################
# STEP 1: Import the raw count matrix 

# searched on the pubmed for scRNAseq data anaylsis, referred to the NCBI Geo dataset GSM5155198, downloaded the data for osteosarcoma.
# If data downloaded in other format for eg. "tar.gz" format, then first covert it into count matrix format using cell ranger pipeline.
# If the GEO sample data is in 10X format with barcodes/genes/matrix, use Read10X(). 
# the name of your files downloaded should be "barcodes.tsv.gz", "features.tsv.gz" and "matrix.mtx.gz" for the Read10X to recognize them. hence rename them first if not.
# barcodes are the cells
# In this step, the count matrix is read and uploaded here as "raw_counts".
#########################################################################
list.files("C:/Users/Acer/Desktop/Veena_Patil_Lab/scRNAseq_analysis_practise/Data/patient_4", all.files = TRUE)

setwd("C:/Users/Acer/Desktop/Veena_Patil_Lab/scRNAseq_analysis_practise/Data/patient_4")

file.rename("GSM5155198_OS_4_barcodes.tsv.gz", "barcodes.tsv.gz")
file.rename("GSM5155198_OS_4_features.tsv.gz", "features.tsv.gz")
file.rename("GSM5155198_OS_4_matrix.mtx.gz",   "matrix.mtx.gz")

# If above steps are already done, we can directly start from here
raw_counts <- Read10X(data.dir = "C:/Users/Acer/Desktop/Veena_Patil_Lab/scRNAseq_analysis_practise/Data/patient_4")

##################################################
# STEP 2: Create a Seurat object
# We create a Seurat object with unfiltered or raw_counts (from count matrix)
# Parameters 'min.cells = 3' is to include gene detected in at least 3 cells
# Parameters 'min.features = 200' is to include cells with at least 200 detected genes
# Parameters filters out the low expressed genes, ambient RNA or empty droplets
# view and range functions are to keep track of the changes in the metadata for every small steps.
##################################################

osteosarcoma <- CreateSeuratObject(counts = raw_counts, project = "Osteosarcoma_scRNAseq", min.cells = 3, min.features = 200)

View(osteosarcoma@meta.data)
range(osteosarcoma$nFeature_RNA)
range(osteosarcoma$nCount_RNA)
#nFeature_RNA: Number of unique genes detected per cell
#nCount_RNA: Total UMI counts per cell

##################################################
# STEP 3: Calculate % mitochondrial reads and add "percent.mt" column to the metadata
# High mitochondrial percentage indicates: Stressed or dying cells, Poor cell viability
# 'percent.mt' represents mitochnodrial gene expression per cell
##################################################

osteosarcoma[["percent.mt"]] <- PercentageFeatureSet(osteosarcoma, pattern = "^MT-")
View(osteosarcoma@meta.data)
##################################################
# STEP 4: Visualize Quality Control (QC) metrics
# violin plot is to visualise the QC metrics across the cells. 
# This futher helps to visualise the outliers and and set filtering thresholds
##################################################

VlnPlot(osteosarcoma, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

FeatureScatter(osteosarcoma, "nFeature_RNA", "percent.mt") 
#Typically shows negative correlation: high-quality cells have moderate gene counts and low mitochondrial percentage.
##################################################
# STEP 5: Filter out low-quality cells
# empty droplets could  have low features and counts. while cell doublets (multiplets) would have higher values for same.
# We keep cells with enough cells, genes and low mitochondrial content 
# Here, I have kept the threshold values exactly same as mentioned in materials and method section of the article
##################################################

osteosarcoma <- subset(osteosarcoma, subset = nFeature_RNA > 500 & nFeature_RNA < 10000 & nCount_RNA < 150000 & percent.mt < 25)

range(osteosarcoma$nFeature_RNA)
range(osteosarcoma$nCount_RNA)
range(osteosarcoma$percent.mt)
#nFeature_RNA > 500: Removes empty droplets/low-quality cells
#nFeature_RNA < 10000: Removes doublets/multiplets
#nCount_RNA < 150000: Removes cells with abnormally high counts
#percent.mt < 25: Removes stressed/dying cells
##################################################
# STEP 6: Normalize the data
# purpose: Corrects for differences in sequencing depth between cells.
# Converts raw counts into normalized expression levels
# take a look at the data to ensure the normalization
##################################################
osteosarcoma <- NormalizeData(osteosarcoma, normalization.method = "LogNormalize", scale.factor = 10000)

osteosarcoma@assays[["RNA"]]$data@x
# Method: Log-normalization. Counts are divided by total counts per cell. Multiplied by scale factor (10,000). Natural log transformation: ln(1 + normalized_count)
# Why log transform?: Stabilizes variance. Makes data more normally distributed. Reduces impact of highly expressed genes

##################################################
# STEP 7: Identify highly variable genes
# Purpose: Identifies genes that show high cell-to-cell variation
# Why important?: These genes drive biological differences between cell types
# Method: "vst" (variance stabilizing transformation). Models mean-variance relationship. Selects top 2000 variable genes
##################################################
osteosarcoma

osteosarcoma <- FindVariableFeatures(osteosarcoma, selection.method = "vst", nfeatures = 2000)

osteosarcoma

#take a look at the variable features using a plot. this visualises mean expression vs variance
VariableFeaturePlot(osteosarcoma)

##################################################
# STEP 8: Data Scaling
# Purpose: Centers and scales gene expression values.
# Why scale?: Prevents highly expressed genes from dominating PCA
# Equal contribution from all genes. Required for dimensionality reduction
# Transformation: For each gene: (expression - mean)/standard deviation
# Note: Scaling all genes is memory-intensive. Often done only on variable genes
##################################################

all_genes <- rownames(osteosarcoma)
osteosarcoma <- ScaleData(osteosarcoma, features = all_genes)

##################################################
# STEP 9: PCA (linear dimensionality reduction)
# Reduces thousands of genes into principal components for better visualisation
# Purpose: Reduces dimensionality while preserving variance.
# Uses only variable genes (computationally efficient). Calculates 50 principal components.
##################################################

osteosarcoma <- RunPCA(osteosarcoma, features = VariableFeatures(osteosarcoma), npcs = 50)

# Shows top genes contributing to each PC.
print(osteosarcoma[["pca"]], dims = 1:5, nfeatures = 5)

# Visualizes gene loadings for PCs 1 and 2
VizDimLoadings(osteosarcoma, dims = 1:2, reduction = "pca")

# examine and view PCA results in a different ways 
# for eg. DimPlot(), VizDimReduction(), and DimHeatmap()
# Dimplot: Projects cells in 2D PCA space.
DimPlot(osteosarcoma, reduction = "pca", dims = c(1,2))

# Elbow plot ranks the PCs based on their standard deviation.
# Purpose: Determines how many PCs to use for downstream analysis.

ElbowPlot(osteosarcoma, ndims = 50)

##################################################
# STEP 10: Finding neighbors and clusters
# The cluster of the cells are formed based on the PCA scores
# values mentioned are followed from the article
# Purpose: Builds a k-nearest neighbor graph.

# Parameters: dims = 1:20: Uses first 20 PCs. k.param = 20: 20 nearest neighbors
# Creates a shared nearest neighbor (SNN) graph for clustering.
##################################################

osteosarcoma <- FindNeighbors(osteosarcoma, dims = 1:20, k.param = 20)
osteosarcoma <- FindClusters(osteosarcoma, resolution = 0.5, algorithm = 1)

##################################################
# STEP 11: non-linear dimensional reduction (UMAP/tSNE)
# Displays cells in 2D, based on their transcriptomes 
##################################################
# values mentioned are followed from the article
# UMAP: Uniform Manifold Approximation and Projection. Preserves both local and global structure better than t-SNE.
# Parameters:
#n.neighbors: Local neighborhood size (balances local/global structure)
#min.dist: Minimum distance between points (controls clustering)


osteosarcoma <- RunUMAP(osteosarcoma, dims = 1:20, n.neighbors = 30, min.dist = 0.3)
DimPlot(osteosarcoma, reduction = "umap", label = TRUE, repel = TRUE)

# t-SNE: t-Distributed Stochastic Neighbor Embedding. Emphasizes local structure over global structure.
osteosarcoma <- RunTSNE(object = osteosarcoma, dims = 1:20)
DimPlot(object = osteosarcoma, reduction = "tsne", label = TRUE)

# for the paper style UMAP
DimPlot(osteosarcoma,reduction = "umap", split.by = "seurat_clusters",
  ncol = 3, pt.size = 0.6) + theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 10))

##################################################
# STEP 12: Find marker genes for clusters
# These help us identify what cell types each cluster represents
# Purpose: Identifies differentially expressed genes for each cluster.
# Parameters:
#only.pos = TRUE: Only positive markers (higher in cluster vs others)
#min.pct = 0.25: Gene expressed in at least 25% of cells in cluster
#logfc.threshold = 0.25: Minimum log2 fold change
#Uses Wilcoxon rank sum test by default.
##################################################

os_markers <- FindAllMarkers(osteosarcoma, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# tried to get the paper style heatmap, but its more like tutorial one 
# Creates a heatmap of top marker genes for each cluster
top10 <- os_markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
DoHeatmap(osteosarcoma, features = top10$gene, group.by = "seurat_clusters", size = 4)

#Name the cell types according to what the article mentions and plot it back on the UMAP/tSNE
Idents(osteosarcoma) <- "seurat_clusters"


# Defines known marker genes for expected cell types in osteosarcoma.
canonical.markers <- list(
  "Inflammatory Myeloid" = c("LYZ", "S100A8", "S100A9", "CTSD"),
  "Macrophages" = c("C1QC", "C1QB", "APOE", "LGALS3"),
  "Neutrophils" = c("FCGR3B", "CSF3R", "MNDA"),
  "T Cells" = c("CD3D", "CD3E", "TRBC1"),
  "NK Cells" = c("NKG7", "GNLY", "KLRD1"),
  "B Cells" = c("CD79A", "MS4A1", "CD37"),
  "Fibroblasts" = c("COL1A1", "COL1A2", "DCN"),
  "Vascular Smooth Muscle Cells" = c("ACTA2", "TAGLN", "MYH11"),
  "Mast Cells" = c("TPSAB1", "TPSB2", "KIT")
)

# Visualizes expression of canonical markers across clusters.
# make a dot plot to ensure the canonical markers expression and confirm the cell types
DotPlot(
  osteosarcoma,
  features = canonical.markers) + RotatedAxis()

#Name the clusters 
#Manually annotates clusters based on marker expression.
#Matches clusters to known cell types from literature.
new.cluster.ids <- c(
  "Macrophages",        # 0
  "NK Cells",           # 1
  "T Cells",            # 2
  "Macrophages",        # 3
  "Fibroblasts",        # 4
  "Neutrophils",        # 5
  "B Cells",            # 6
  "Inflammatory Myeloid", # 7
  "VSMC",               # 8
  "Mast Cells",         # 9
  "Inflammatory Myeloid", # 10
  "T Cells",            # 11
  "Fibroblasts",        # 12
  "Macrophages"         # 13
)

names(new.cluster.ids) <- levels(osteosarcoma)
osteosarcoma <- RenameIdents(osteosarcoma, new.cluster.ids)

DimPlot(osteosarcoma, reduction = "umap", label = TRUE, repel = TRUE)
DimPlot(osteosarcoma, reduction = "tsne", label = TRUE, repel = TRUE)

#violin plot for NFATC1 expression across the cells
VlnPlot(osteosarcoma, features = "NFATC1", pt.size = 0.1)
FeaturePlot(osteosarcoma, features = "NFATC1")

#for NFATc1 expression in Vascular Smooth Muscle Cells (VMSCs)
vsmc <- subset(osteosarcoma, idents = "VSMC")
FeaturePlot(vsmc, features = "NFATC1", reduction = "umap", pt.size = 1.2)


##################################################
# STEP 14: Save formatted Seurat object for future use
##################################################

saveRDS(osteosarcoma, file = "osteosarcoma_scRNA_final.rds")
