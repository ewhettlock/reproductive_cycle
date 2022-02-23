#The endometrial scRNAseq data set (available at www.reproductivecellatlas.org/) was avaliable in the python format so initially converted into an R Seurat object

#Loading libraries for conversion from python dataset to R dataset using conda
library(rhdf5)
library(reticulate)

# Create a new environment 
conda_create("endor")

# install anndata
conda_install("endor", "anndata")

#making h5 into seurat object for use in R
library(rio)
endo <- Convert("*/endometrium_all.h5ad", dest = "h5seurat", overwrite = TRUE)
endo <- LoadH5Seurat(endo)

#Confirming "endo" is now a seurat object and saving the .rds file
endo
saveRDS(endo, file = "*/endo_original.rds")

#Loading .rds file
endo <- readRDS(file = "*/endo_original.rds")

#Loading required libraries for seurat work
library(Seurat)
library(dplyr)
library(patchwork)


#Onto to analysis - checking how current data looks
head(endo@meta.data, 5)
VlnPlot(endo, features = c("n_genes", "log2p1_count", "percent_mito"), ncol = 3)

#To look at data by itself can either use SCTransform or older method (however if just using to combine datasets then skip these)

#Method 1
SCTransform(endo)

#Method 2
endo <- NormalizeData(endo)
endo <- FindVariableFeatures(endo)
endo <- ScaleData(endo)

#Both methods join here to run PCA analysis and visualise on a UMAP
endo <- RunPCA(endo, features = VariableFeatures(object = endo))
ElbowPlot(endo)
endo <- FindNeighbors(endo, dims = 1:18)
endo <- FindClusters(endo, resolution = 0.5)
endo <- RunUMAP(endo, dims = 1:18)
DimPlot(endo, reduction = "umap")
DimPlot(endo, reduction = "umap", group.by = "Cell.type")

#Can also look at features of the data e.g.
FeaturePlot(endo, features = c("CD3E", "NKG7"))
VlnPlot(endo, features = c("CD3E", "NKG7"))

#However since this dataset has been analysed and already has meta data which assigns cell.type, the above analysis is not required
# I was only interested in immune cells in the endomterium so I selected cells of the lymphoid and myeloid type, and removed cells from the myometrium
endo_immune <- subset(endo, subset = Cell.type == "Lymphoid" | Cell.type == "Myeloid")
endo_immune_nomyo <- subset(endo_immune, subset = Location != "myometrium")
saveRDS(endo_immune_nomyo, file = "*/endo_immune_nomyo.rds")



