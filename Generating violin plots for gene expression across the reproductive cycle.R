#Loading libraries
library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)

#setting working directory
setwd("*")

#Loading integrated data
all <- readRDS("FourData_paper.rds")

#Subsetting to individual datasets for the downsampling
non_preg <- subset(all, Dataset == "non preg")
first <- subset(all, Dataset == "first")
second <- subset(all, Dataset == "second")
third <- subset(all, Dataset == "third")

#Whole object subsetting so each dataset has the same number of cells for visualisation by dataset in violin plots
n <- sample(colnames(non_preg), size = 2000, replace = F)
f <- sample(colnames(first), size = 2000, replace = F)
s <- sample(colnames(second), size = 2000, replace = F)
t <- sample(colnames(third), size = 2000, replace = F)
all <- all[, c(n, f, s, t)]

#Subsetting to NK and changing labelling
NK <- subset(all, subset = seurat_clusters == 3 | seurat_clusters == 1 | seurat_clusters == 5)
new.cluster.ids <- c("dNK2", "dNK1", "dNK3")
names(new.cluster.ids) <- levels(NK)
NK <- RenameIdents(NK, new.cluster.ids)
DimPlot(NK, reduction = "umap", label = TRUE)

#Adding the cell identity to the metadata
dNK <- Idents(object = NK)
names(dNK) <- colnames(NK)
NK <- AddMetaData(object = NK, metadata = dNK, col.name = "dNK_type")
NK$dNK_type <- factor(NK$dNK_type, levels = c("dNK1", "dNK2", "dNK3"))

#Switching back to RNA to make the violin plots
DefaultAssay(NK) <- "RNA"

#Change directory to new location for saving violin plots
setwd("*/VlnPlots")

#Selecting the chosen colours
purple <- rgb(red = (118/255), green = (6/255), blue = (154/255))
orange <- rgb(red = (225/255), green = (96/255), blue = (0/255))
green <- rgb(red = (0/255), green = (128/255), blue = (0/255))

#Removing old datasets to keep environment tidy
rm(all)
rm(non_preg)
rm(first)
rm(second)
rm(third)

#Function for producing violin plots for expression of gene "x". 3 plots made for dNK1, 2 and 3.
#Each plot has the same y limit for easier comparison
#Each plot has a different colour dNK1 = purple, dNK2 = orange, dNK3 = green

Plot <- function(x) {
  
  limit = FetchData(NK, vars = x) %>% max() %>% ceiling()
  
  stats <- FetchData(NK, vars = c("Time", x , "dNK_type"))
  stats1 <- stats %>% subset(dNK_type == "dNK1")
  
  stats1 %>%
    ggplot(aes(x = Time, y = stats1[,2], color = dNK_type, fill = dNK_type)) +
    geom_violin() +
    ylim(0, limit) +
    scale_colour_manual(values = purple) +
    scale_fill_manual(values = purple) +
    labs(title = "", y = "Expression level") +
    theme(legend.position = 'none',
          text = element_text(size = 8),
          axis.text.x = element_text(size = 8, angle=50, hjust=1),
          axis.text.y = element_text(size = 8),
          panel.background = element_rect(fill = "NA"),
          axis.line = element_line(colour = "black"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 8))
  
  ggsave(paste0(x, "_dNK1", ".pdf"), width = 5, height = 5, units = "cm")
  
  stats2 <- stats %>% subset(dNK_type == "dNK2")
  
  stats2 %>%
    ggplot(aes(x = Time, y = stats2[,2], color = dNK_type, fill = dNK_type)) +
    geom_violin() +
    ylim(0, limit) +
  scale_colour_manual(values = orange) +
    scale_fill_manual(values = orange) +
    labs(title = "", y = "Expression level") +
    theme(legend.position = 'none',
          text = element_text(size = 8),
          axis.text.x = element_text(size = 8, angle=50, hjust=1),
          axis.text.y = element_text(size = 8),
          panel.background = element_rect(fill = "NA"),
          axis.line = element_line(colour = "black"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 8))
  
  ggsave(paste0(x, "_dNK2", ".pdf"), width = 5, height = 5, units = "cm")
  
  stats3 <- stats %>% subset(dNK_type == "dNK3")
  
  stats3 %>%
    ggplot(aes(x = Time, y = stats3[,2], color = dNK_type, fill = dNK_type)) +
    geom_violin() +
    ylim(0, limit) +
  scale_colour_manual(values = green) +
    scale_fill_manual(values = green) +
    labs(title = "", y = "Expression level") +
    theme(legend.position = 'none',
          text = element_text(size = 8),
          axis.text.x = element_text(size = 8, angle=50, hjust=1),
          axis.text.y = element_text(size = 8),
          panel.background = element_rect(fill = "NA"),
          axis.line = element_line(colour = "black"),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 8))
  
  
  ggsave(paste0(x, "_dNK3", ".pdf"), width = 5, height = 5, units = "cm")
}

#Example for one gene
Plot("LILRB1")

#Used a gene vector and sapply to plot multiple genes at once
Genes <- c("CSF2", "GNLY", "GZMA", "GZMB", "GZMH", "GZMK", "GZMM", "IFNG", "KIR2DL1", "KIR2DL3", "KIR2DL4", "KLRD1", "LAMP1", "LILRB1", "PRF1", "TNF", "XCL1", "CXCL8", "VEGFA", "VEGFB", "VEGFC")
sapply(X = Genes, FUN = Plot)

#This function was also used to generate the labouring violin plots

