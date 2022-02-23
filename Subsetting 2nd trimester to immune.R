#The scRNAseq dataset (dbGaP phs001886.v1.p1, reanalyzed with permission of the NIH, project ID 145 26528) contained one accreta sample
#This data was given in .sra files. This was converted to fastq files and run through cell ranger using h19 as a reference.

#Loading libraries
library(dplyr)
library(Seurat)
library(patchwork)

#Loading the data from cell ranger. Excludes features that appear in less than 10 cells. Excludes cells that have less than 200 features
second.data <- Read10X(data.dir = "*/filtered_feature_bc_matrix")
second <- CreateSeuratObject(counts = second.data, min.cells = 10, min.features = 200)

#Confirming what the seurat object looks like and looking at the meta.data
second
head(second@meta.data)

#Adding mitochondrial data to the metadata
second[["percent.mt"]] <- PercentageFeatureSet(second, pattern = "^MT-")

#Exclude cells where over 10% of the genes are mitochondrial
second <- subset(second, subset = nFeature_RNA > 200 & percent.mt < 10)

#SCtransform 
second <-SCTransform(second)

#PCA
second <- RunPCA(second, features = VariableFeatures(object = second))
DimPlot(second, reduction = "pca")

#Use elbow plot to etermine dimensions used in the following algorithms
ElbowPlot(second)

#Clustering
second <- FindNeighbors(second, dims = 1:20)
second <- FindClusters(second, resolution = 0.5)

#Visualising data using UMAP
second <- RunUMAP(second, dims = 1:20)
DimPlot(second, reduction = "umap")
DimPlot(second, reduction = "umap", label = TRUE, pt.size = 0.5)

#Finding biomarkers for each cluster and then generating a top 10 list
cluster.markers <- FindAllMarkers(second, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- cluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

#Imaging individual gene expression on umap
FeaturePlot(second, features = "PTPRC")
FeaturePlot(second, features = c("GNLY", "CD3E", "CD14", "PAGE4", "FSTL3", "RGS5", "DKK1", "DLK1"))


#Used biomarkers to identify which clusters were immune and which were non-immune. Subsetting to just immune clusters
Pique_Immune <- subset(piqueB, subset = seurat_clusters == 0 | seurat_clusters ==  4 | seurat_clusters == 6 | seurat_clusters == 7 | 
                         seurat_clusters == 9 | seurat_clusters == 10 | seurat_clusters == 13 | seurat_clusters == 15 | 
                         seurat_clusters == 16 | seurat_clusters == 17 | seurat_clusters == 19 | seurat_clusters == 21 | seurat_clusters == 22)

saveRDS(Pique_Immune, file = "*/2nd_immune.rds")


