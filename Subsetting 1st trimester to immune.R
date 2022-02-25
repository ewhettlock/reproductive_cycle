#First trimester dataset (available at Array Express E-MTAB-6701) was converted from a .txt file to an R Seurat object.
#Then subsetted to immune cells using the "annotation" column of the metadata
# This code was written with the help of this tutorial https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

#Loading libraries
library(dplyr)
library(Seurat)
library(patchwork)

#Loading data in the form of a .txt file
VT <- read.delim("raw_data.txt")

#Appears that first column has the gene  names. Need to make the gene names the rownames instead, in order to convert to Seurat object
dim(VT)
head(rownames(VT))

#Select gene names from first column into a vector "G"
G <- VT$Gene

#Editing to remove ensemble numbers to be consistent with gene names from other scRNAseq data sets
H <- sapply(X = strsplit(G, split = "_EN"), FUN = "[", 1)

#The Gene name column has some dulpicates. Adding -1 to one of each duplicate so all gene names are unique as required
n_occur <- data.frame(table(H))
n_occur[n_occur$Freq > 1,]

H[match("ATRIP", H)] <- "ATRIP-1"
H[match("BLOC1S5", H)] <- "BLOC1S5-1"
H[match("CCDC7", H)] <- "CCDC7-1" 
H[match("CFAP99", H)] <- "CFAP99-1" 
H[match("CYB561D2", H)] <- "CYB561D2-1" 
H[match("IGHV2-70", H)] <- "IGHV2-70-1"
H[match("LINC00864", H)] <- "LINC00864-1"
H[match("LINC01297", H)] <- "LINC01297-1" 
H[match("LINC01422", H)] <- "LINC01422-1"
H[match("LINC01481", H)] <- "LINC01481-1" 
H[match("MATR3", H)] <- "MATR3-1"
H[match("PGM5-AS1", H)] <- "PGM5-AS1-1" 
H[match("PRICKLE4", H)] <- "PRICKLE4-1"
H[match("RAET1E-AS1", H)] <- "RAET1E-AS1-1" 
H[match("RGS5", H)] <- "RGS5-1" 
H[match("SERPINA3", H)] <- "SERPINA3-1" 
H[match("SPATA13", H)] <- "SPATA13-1" 
H[match("TBC1D26", H)] <- "TBC1D26-1" 
H[match("TIMM10B", H)] <- "TIMM10B-1" 
H[match("TMEM256-PLSCR3", H)] <- "TMEM256-PLSCR3-1" 

n_occur <- data.frame(table(H))
n_occur[n_occur$Freq > 1,]


#Make Gene names the rownames
row.names(VT) <- H

#Remove the gene column from the data frame
VT <- subset(VT, select = -Gene)

#Checking the above has worked
dim(VT)
head(rownames(VT))

#Create Seurat Object
VT <- CreateSeuratObject(VT, min.cells = 3, min.features = 500)

#Adding meta data
meta <- read.delim("meta_10x.txt")
VT <- AddMetaData(VT, meta)

#Checking meta data
head(VT@meta.data)
tail(VT@meta.data)

#SCTransform, PCA and UMAP analysis
VT <- SCTransform(VT)
VT <- RunPCA(VT, features = VariableFeatures(object = VT))
VT <- FindNeighbors(VT, dims = 1:20)
VT <- FindClusters(VT, resolution = 0.5)
VT <- RunUMAP(VT, dims = 1:20)

saveRDS(VT, file = "VT_Meta_SCT_UMAP.rds")

#Subsetting to just decidual cells (removing blood and placental cells)
VT <- subset(VT, subset = location == "Decidua")

#Selecting immune cells using the "annotation" column of the meta data
VT_Immune <- subset(VT, subset = annotation == "Granulocytes" | annotation == "Plasma" | annotation == "DC1" | annotation == "DC2" | annotation == "dM1" | annotation == "dM2" | annotation == "dM3" | annotation == "dNK p" | annotation == "dNK1" | annotation == "dNK2" | annotation == "dNK3" | annotation == "HB" | annotation == "ILC3" | annotation == "MO" | annotation == "NK CD16-" | annotation == "NK CD16+" | annotation == "Tcells")

#Saving the .rds file
saveRDS(VT_Immune, file = "VT_Immune_Meta_Ann.rds")
