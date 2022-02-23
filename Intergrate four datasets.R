#Loading libraries
library(dplyr)
library(Seurat)
library(patchwork)

#Loading the four individual datasets
second <- readRDS(file = "2nd_immune.rds")
endo <- readRDS(file = "endo_immune_nomyo.rds")
third <- readRDS(file = "pique_immuneUMAP.rds")
first <- readRDS(file= "VT_Immune_Meta_Ann.rds")

#Confirming the default assay for all is "RNA"
DefaultAssay(second) <- "RNA"
DefaultAssay(endo) <-	"RNA"
DefaultAssay(third) <-	"RNA"
DefaultAssay(first) <-	"RNA"

#Removing placental and labouring data from the third trimester dataset 
third <- subset(third, subset = Tissue != "Placental villi")
third <- subset(third, subset = Condition == "TNL")

#Removing a third trimester particpiant that did not overlap with the other cells
first <- subset(first, subset = Fetus != "D12")

#Merging into one object
all <- merge(endo, y = c(first, second, third), add.cell.ids = c("non preg", "first", "second", "third"))

#Adding dataset to meta data
all[["Dataset"]] <- sapply(X = strsplit(colnames(all), split = "_"), FUN = "[", 1)

#split by dataset
all.list <- SplitObject(all, split.by = "Dataset")

#Normalise and identify variable features for each set independently 
all.list <- lapply(X=all.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  y <- x@assays$RNA@var.features[!grepl('^MT|^RP', x@assays$RNA@var.features)]
  VariableFeatures(x) <- y
  x
})

#select feature that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = all.list)

#Performing intergration
immune.anchors <- FindIntegrationAnchors(object.list = all.list, anchor.features = features)
immune.combined <- IntegrateData(anchorset = immune.anchors)

#Integrated analysis
DefaultAssay(immune.combined) <- "integrated"

#Standard workflow
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:30)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)

#Saving combinded file
saveRDS(immune.combined, file = "FourData.rds")

