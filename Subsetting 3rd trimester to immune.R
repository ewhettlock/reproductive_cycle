#The scRNAseq dataset (dbGaP phs001886.v1.p1, reanalyzed with permission of the NIH, project ID 145 26528) contained nine third trimester partipants
#This data was given in .sra files. This was converted to fastq files and run through cell ranger using h19 as a reference.

#Loading libraries
library(dplyr)
library(Seurat)
library(patchwork)

#A seurat object was created using min.cells = 10 and min.feature = 200 when the meta data was added
pique <- readRDS(file = "*/pique_with_meta.rds")

##Adding mitochondrial data to the metadata
pique[["percent.mt"]] <- PercentageFeatureSet(pique, pattern = "^MT-")

#Exclude cells where over 10% of the genes are mitochondrial
pique <- subset(pique, subset = nFeature_RNA > 200 & percent.mt < 10)

#SCTransform
pique <- SCTransform(pique)

#Removing mitochondrial and ribosomal genes from the variable features lists
VariableFeatures(pique) <- pique@assays$SCT@var.features[!grepl('^MT|^RP', pique@assays$SCT@var.features)]

#PCA analysis and UMAP
pique <- RunPCA(pique, features = VariableFeatures(object = pique))
pique <- FindNeighbors(pique, dims = 1:20)
pique <- FindClusters(pique, resolution = 0.5)
pique <- RunUMAP(pique, dims = 1:20)

#Finding Cluster Markers to help identify clusters
cluster.markers <- FindAllMarkers(pique, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- cluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

#subsetting to just immune cells based on the cluster identities from the cluster markers
pique_immune <- subset(pique, subset = seurat_clusters == 0 | seurat_clusters ==  3 | seurat_clusters == 4 | seurat_clusters == 5 | 
                  seurat_clusters == 7 | seurat_clusters == 8 | seurat_clusters == 13 | seurat_clusters == 16)


saveRDS(pique_Immune, file = "*/pique_immuneUMAP.rds")

