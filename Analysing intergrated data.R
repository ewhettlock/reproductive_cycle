#Loading libraries
library(Seurat)
library(patchwork)
library(dplyr)

#setting working directory
setwd("*")

#Loading integrated data
all <- readRDS("Fourdata.rds")

#Finding biomarkers for each cluster and then generating a top 10 list
cluster.markers <- FindAllMarkers(all, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- cluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

#Changing resolution and confirming by UMAP that resolution is okay
all <- FindClusters(all, resolution = 0.48)
DimPlot(all, reduction = "umap", label = TRUE)

#Checking meta data
head(all@meta.data)
tail(all@meta.data)

#Changing order of Datasets to order in reproductive cycle
all$Dataset <- factor(all$Dataset, levels = c("non preg", "first", "second", "third"))

#Looking at the number of cells from each dataset
cell_numbers <- table(all$Dataset)

#Adding stages for different datasets into one consistent column in the 
Tissue <- all$Tissue
Binary_Stage <- all$Binary.Stage
Dataset <- all$Dataset
Meta <- data.frame(Tissue, Binary_Stage, Dataset)

Add2 <- function(x, y, z) {
  if (z == "first") {
    print("1T")
  } else if (z == "second") {
    print("2T")
  } else if (z == "non preg") {
    print(y)
  } else if (z == "third") {
    print(x)
  } else print("Help")
}

Meta$Time <- mapply(x = Meta$Tissue, y = Meta$Binary_Stage, z = Meta$Dataset, FUN = Add2)

Meta$Time <- gsub("Basal Plate", "3T (DB)", Meta$Time)
Meta$Time <- gsub("Chorioamniotic membranes", "3T (DP)", Meta$Time)

Meta <- subset(Meta, select = "Time")

all <- AddMetaData(all, Meta)

#Changing order of new variable "Time" to order in the reproductive cycle
all$Time <- factor(all$Time, levels = c("Proliferative", "Secretory", "1T", "2T", "3T (DB)", "3T (DP)"))

saveRDS(all, file = "Fourdata_paper.rds")

# # #
# Adding consistent sample IDs across datasets
Third <- all$SampleID
First <- all$Fetus
NonPreg <- all$DonorID
Meta <- data.frame(NonPreg, First, Third)

Add2 <- function(x, y, z) {
  if (is.na(x) == "TRUE" & is.na(y) == "TRUE" & is.na(z) == "TRUE") {
    print("2T")
  } else if (is.na(x) == "TRUE" & is.na(y) == "TRUE") {
    print(z)
  } else if (is.na(x) == "TRUE" & is.na(z) == "TRUE") {
    print(y)
  } else print(x)
}

Meta$Sample_Tag <- mapply(x = Meta$NonPreg, y = Meta$First, z = Meta$Third, FUN = Add2)

Meta <- subset(Meta, select = "Sample_Tag")


# #
# #Changing sample IDs to be consistent
#
Meta$Sample_Tag <- gsub("D6", "First_D6", Meta$Sample_Tag)
Meta$Sample_Tag <- gsub("D7", "First_D7", Meta$Sample_Tag)
Meta$Sample_Tag <- gsub("D8", "First_D8", Meta$Sample_Tag)
Meta$Sample_Tag <- gsub("D9", "First_D9", Meta$Sample_Tag)
Meta$Sample_Tag <- gsub("D10", "First_D10", Meta$Sample_Tag)
Meta$Sample_Tag <- gsub("A13", "NonPreg_A13", Meta$Sample_Tag)
Meta$Sample_Tag <- gsub("A30", "NonPreg_A30", Meta$Sample_Tag)
Meta$Sample_Tag <- gsub("E1", "NonPreg_E1", Meta$Sample_Tag)
Meta$Sample_Tag <- gsub("E2", "NonPreg_E2", Meta$Sample_Tag)
Meta$Sample_Tag <- gsub("E3", "NonPreg_E3", Meta$Sample_Tag)
Meta$Sample_Tag <- gsub("SAMN15049042", "NonPreg_SAM42", Meta$Sample_Tag)
Meta$Sample_Tag <- gsub("SAMN15049043", "NonPreg_SAM43", Meta$Sample_Tag)
Meta$Sample_Tag <- gsub("SAMN15049044", "NonPreg_SAM44", Meta$Sample_Tag)
Meta$Sample_Tag <- gsub("SAMN15049045", "NonPreg_SAM45", Meta$Sample_Tag)
Meta$Sample_Tag <- gsub("SAMN15049046", "NonPreg_SAM46", Meta$Sample_Tag)
Meta$Sample_Tag <- gsub("SAMN15049047", "NonPreg_SAM47", Meta$Sample_Tag)
Meta$Sample_Tag <- gsub("SAMN15049048", "NonPreg_SAM48", Meta$Sample_Tag)
Meta$Sample_Tag <- gsub("SAMN15049049", "NonPreg_SAM49", Meta$Sample_Tag)
Meta$Sample_Tag <- gsub("SAMN15049050", "NonPreg_SAM50", Meta$Sample_Tag)
Meta$Sample_Tag <- gsub("SAMN15049051", "NonPreg_SAM51", Meta$Sample_Tag)

#To check the above as worked
table(Meta$Sample_Tag)

# Adding new column to the meta data 
all <- AddMetaData(all, Meta)

saveRDS(all, file = "Fourdata_paper.rds")

#Visulalising by four datasets
DimPlot(all, reduction = "umap", split.by = "Dataset")

#Subsetting to indivdual datasets
non_preg <- subset(all, subset = Dataset == "non preg")
first <- subset(all, subset = Dataset == "first")
second <- subset(all, subset = Dataset == "second")
third <- subset(all, subset = Dataset == "third")

#Whole object subsetting so each dataset has the same number of cells for visualisation by dataset
n <- sample(colnames(non_preg), size = 2000, replace = F)
f <- sample(colnames(first), size = 2000, replace = F)
s <- sample(colnames(second), size = 2000, replace = F)
t <- sample(colnames(third), size = 2000, replace = F)
all_2000 <- all[, c(n, f, s, t)]

#Visualising cells by dataset, with each dataset having an equal number of cells
DimPlot(all_2000, reduction = "umap", split.by = "Dataset", pt.size = 0.8)

#Finding biomarkers for each cluster and then generating a top 10 list
cluster.markers <- FindAllMarkers(all, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- cluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

#Need to remove 10 and 15 and 19 and 22 as they do not appear to be immune
all <- subset(all, subset = seurat_clusters != 10)
all <- subset(all, subset = seurat_clusters != 15)
all <- subset(all, subset = seurat_clusters != 19)
all <- subset(all, subset = seurat_clusters != 22)


# Rerunning analysis due to removal of non immune clusters
DefaultAssay(all) <- "integrated"
all <- RunPCA(all, verbose = FALSE)
all <- FindNeighbors(all, reduction = "pca", dims = 1:30)
all <- FindClusters(all, resolution = 0.33)
all <- RunUMAP(all, reduction = "pca", dims = 1:30)
DimPlot(all, reduction = "umap", label = TRUE)
saveRDS(all, file = "Fourdata_paper.rds")

#Switching back to RNA for looking at gene expression
DefaultAssay(all) <- "RNA"

#Finding biomarkers for each cluster in new set and then generating a top 10 list
cluster.markers <- FindAllMarkers(all, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- cluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

#Feature Plots
FeaturePlot(all, features = c("NKG7", "CD3E"))

#Violin Plots
VlnPlot(all, features = c("IL2RA"))

#Table
table(all$Dataset, all$seurat_clusters)
table(all$annotation, all$seurat_clusters)




