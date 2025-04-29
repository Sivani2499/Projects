# standard workflow steps to analyze single cell RNA-Seq data
# data: 17k_Ovarian_Cancer_scFFPE_count_filtered_feature_bc_matrix.h5
# data source:https://www.10xgenomics.com/datasets/17k-human-ovarian-cancer-scFFPE 

#setwd("/Users/shivaniravindran/Library/CloudStorage/Box-Box/R_single_cell_workflows/")


# load libraries
library(Seurat)
library(tidyverse)

# Load the NSCLC dataset
#OC refers to ovarian cancer
OC.sparse.m <- Read10X_h5(filename = '/Users/shivaniravindran/Library/CloudStorage/Box-Box/R_single_cell_workflows/17k_Ovarian_Cancer_scFFPE_count_filtered_feature_bc_matrix.h5')
str(OC.sparse.m)
#For S4 object you can check for the slot names by the following function:
slotNames(OC.sparse.m)
#You can access the slot by using the @ operator instead of the $ operator
counts <-  OC.sparse.m@Dimnames


#here since the oc.sparse.m already contains the gene expression data ( the counts) we are using it directly in the blwo step
# Initialize the Seurat object with the raw (non-normalized data).
oc.seurat.obj <- CreateSeuratObject(counts = OC.sparse.m, project = "ovarian_cancer", min.cells = 3, min.features = 200)
str(oc.seurat.obj)
oc.seurat.obj
# 17354 features across 16956 samples


# 1. QC -------
View(oc.seurat.obj@meta.data)
# % Mitochondrail(MT) reads
oc.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(oc.seurat.obj, pattern = "^MT-")
View(oc.seurat.obj@meta.data)

VlnPlot(oc.seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(oc.seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')

# 2. Filtering -----------------
oc.seurat.obj <- subset(oc.seurat.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 1500 & 
                             percent.mt < 5)

# 3. Normalize data ----------
#oc.seurat.obj <- NormalizeData(oc.seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
# OR
oc.seurat.obj <- NormalizeData(oc.seurat.obj)
str(oc.seurat.obj)


# 4. Identify highly variable features --------------
oc.seurat.obj <- FindVariableFeatures(oc.seurat.obj, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(oc.seurat.obj), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(oc.seurat.obj)
LabelPoints(plot = plot1, points = top10, repel = TRUE)


# 5. Scaling -------------
all.genes <- rownames(oc.seurat.obj)
oc.seurat.obj <- ScaleData(oc.seurat.obj, features = all.genes)

str(oc.seurat.obj)

# 6. Perform Linear dimensionality reduction --------------
oc.seurat.obj <- RunPCA(oc.seurat.obj, features = VariableFeatures(object = oc.seurat.obj))

# visualize PCA results
print(oc.seurat.obj[["pca"]], dims = 1:5, nfeatures = 5)
DimHeatmap(oc.seurat.obj, dims = 1, cells = 500, balanced = TRUE)


# determine dimensionality of the data
ElbowPlot(oc.seurat.obj)


# 7. Clustering ------------
oc.seurat.obj <- FindNeighbors(oc.seurat.obj, dims = 1:15)

# understanding resolution
oc.seurat.obj <- FindClusters(oc.seurat.obj, resolution = c(0.1,0.3, 0.5, 0.7, 1))
View(oc.seurat.obj@meta.data)

DimPlot(oc.seurat.obj, group.by = "RNA_snn_res.1", label = TRUE)

# setting identity of clusters
Idents(oc.seurat.obj)
Idents(oc.seurat.obj) <- "RNA_snn_res.0.1"
Idents(oc.seurat.obj)

# non-linear dimensionality reduction --------------
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
reticulate::py_install(packages ='umap-learn')
oc.seurat.obj <- RunUMAP(oc.seurat.obj, dims = 1:15)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(oc.seurat.obj, reduction = "umap")
