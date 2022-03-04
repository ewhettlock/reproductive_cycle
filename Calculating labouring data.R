#Using the third trimester data set in isolation doesn't cluster the NK cells in to 1,2 and 3.
#However the cells do cluster like that when intergrated with other cells from the reproductive cycle.
#Therefore to see how the frequency of dNK1-3 varies in labour I used the integrated set and subset it to the third trimester dataset

#Loading libraries
library(dplyr)
library(Seurat)
library(patchwork)

#Loading merged dataset 
setwd("*")
all <- readRDS("FourData_paper.rds")

# Selecting just term data
Third <- subset(all, subset = Dataset == "third")

#Visualising term data
DimPlot(Third, reduction = "umap", label = TRUE, pt.size = 0.5)

#Fetching data
#Cluster 5 = dNK1, 1 = dNK2, 2 = dNK3 and 7 = dNKp
stats <- FetchData(Third, vars = c("Condition", "Tissue", "SampleID", "seurat_clusters", "GNLY"))

cell <- c(rep(1, times = 14976))
stats <- mutate(stats, cell)

DB <- stats %>% filter(Tissue == "Basal Plate")
DP <- stats %>% filter(Tissue == "Chorioamniotic membranes")

#Generating stats for the frequency of dNK in the decidua basalis

Stats_summary <- DB %>% 
  filter(seurat_clusters != 3, seurat_clusters != 6, seurat_clusters != 8, seurat_clusters != 12, seurat_clusters != 13, seurat_clusters != 14, seurat_clusters != 15, seurat_clusters != 18) %>%
  group_by(SampleID, Condition) %>%
  arrange(SampleID) %>%
  summarize(Total_lymph = sum(cell))

dNK1 <- DB %>% 
  group_by(SampleID, seurat_clusters) %>%
  summarize(N = n()) %>%
  ungroup() %>%
  complete(SampleID, seurat_clusters, fill = list(N=0)) %>%
  filter(seurat_clusters == 5) %>%
  arrange(SampleID) %>%
  select(N)

Stats_summary$dNK1 <- dNK1

dNK2 <- DB %>% 
  group_by(SampleID, seurat_clusters) %>%
  summarize(N = n()) %>%
  ungroup() %>%
  complete(SampleID, seurat_clusters, fill = list(N=0)) %>%
  filter(seurat_clusters == 1) %>%
  arrange(SampleID) %>%
  select(N)

Stats_summary$dNK2 <- dNK2

dNK3 <- DB %>% 
  group_by(SampleID, seurat_clusters) %>%
  summarize(N = n()) %>%
  ungroup() %>%
  complete(SampleID, seurat_clusters, fill = list(N=0)) %>%
  filter(seurat_clusters == 2) %>%
  arrange(SampleID) %>%
  select(N)

Stats_summary$dNK3 <- dNK3

dNKp <- DB %>% 
  filter(GNLY >= 1) %>%
  group_by(SampleID, seurat_clusters) %>%
  summarize(N = n()) %>%
  ungroup() %>%
  complete(SampleID, seurat_clusters, fill = list(N=0)) %>%
  filter(seurat_clusters == 7) %>%
  arrange(SampleID) %>%
  select(N)

Stats_summary$dNKp <- dNKp

Total_immune <- DB %>% 
  group_by(SampleID, Condition) %>%
  arrange(SampleID) %>%
  summarize(N = sum(cell)) %>%
  ungroup() %>%
  complete(SampleID, fill = list(N=0)) %>%
  arrange(SampleID) %>%
  select(N)


Stats_summary$Total_Immune <- Total_immune

#Using mutate to calculate the percentages

Stats_summary <- Stats_summary %>% mutate(Total_NK = sum(dNK1, dNK2, dNK3, dNKp))
Stats_summary <- Stats_summary %>% mutate(Percent_NK_of_Lymph = (Total_NK/Total_lymph)*100)
Stats_summary <- Stats_summary %>% mutate(Percent_dNK1_of_Lymph = (dNK1/Total_lymph)*100)
Stats_summary <- Stats_summary %>% mutate(Percent_dNK2_of_Lymph = (dNK2/Total_lymph)*100)
Stats_summary <- Stats_summary %>% mutate(Percent_dNK3_of_Lymph = (dNK3/Total_lymph)*100)
Stats_summary <- Stats_summary %>% mutate(Percent_dNK1_of_NK = (dNK1/Total_NK)*100)
Stats_summary <- Stats_summary %>% mutate(Percent_dNK2_of_NK = (dNK2/Total_NK)*100)
Stats_summary <- Stats_summary %>% mutate(Percent_dNK3_of_NK = (dNK3/Total_NK)*100)

#and doing the exact same for the decidua parietalis 

Stats_summaryP <- DP %>% 
  filter(seurat_clusters != 3, seurat_clusters != 6, seurat_clusters != 8, seurat_clusters != 12, seurat_clusters != 13, seurat_clusters != 14, seurat_clusters != 15, seurat_clusters != 18) %>%
  group_by(SampleID, Condition) %>%
  arrange(SampleID) %>%
  summarize(Total_lymph = sum(cell))

dNK1P <- DP %>% 
  group_by(SampleID, seurat_clusters) %>%
  summarize(N = n()) %>%
  ungroup() %>%
  complete(SampleID, seurat_clusters, fill = list(N=0)) %>%
  filter(seurat_clusters == 5) %>%
  arrange(SampleID) %>%
  select(N)

Stats_summaryP$dNK1 <- dNK1P

dNK2P <- DP %>% 
  group_by(SampleID, seurat_clusters) %>%
  summarize(N = n()) %>%
  ungroup() %>%
  complete(SampleID, seurat_clusters, fill = list(N=0)) %>%
  filter(seurat_clusters == 1) %>%
  arrange(SampleID) %>%
  select(N)

Stats_summaryP$dNK2 <- dNK2P

dNK3P <- DP %>% 
  group_by(SampleID, seurat_clusters) %>%
  summarize(N = n()) %>%
  ungroup() %>%
  complete(SampleID, seurat_clusters, fill = list(N=0)) %>%
  filter(seurat_clusters == 2) %>%
  arrange(SampleID) %>%
  select(N)

Stats_summaryP$dNK3 <- dNK3P

dNKpP <- DP %>% 
  filter(GNLY >= 1) %>%
  group_by(SampleID, seurat_clusters) %>%
  summarize(N = n()) %>%
  ungroup() %>%
  complete(SampleID, seurat_clusters, fill = list(N=0)) %>%
  filter(seurat_clusters == 7) %>%
  arrange(SampleID) %>%
  select(N)

Stats_summaryP$dNKp <- dNKpP

Total_immuneP <- DP %>% 
  group_by(SampleID, Condition) %>%
  arrange(SampleID) %>%
  summarize(N = sum(cell)) %>%
  ungroup() %>%
  complete(SampleID, fill = list(N=0)) %>%
  arrange(SampleID) %>%
  select(N)


Stats_summaryP$Total_Immune <- Total_immuneP

#Calculating the percentages using mutate

Stats_summaryP <- Stats_summaryP %>% mutate(Total_NK = sum(dNK1, dNK2, dNK3, dNKp))
Stats_summaryP <- Stats_summaryP %>% mutate(Percent_NK_of_Lymph = (Total_NK/Total_lymph)*100)
Stats_summaryP <- Stats_summaryP %>% mutate(Percent_dNK1_of_Lymph = (dNK1/Total_lymph)*100)
Stats_summaryP <- Stats_summaryP %>% mutate(Percent_dNK2_of_Lymph = (dNK2/Total_lymph)*100)
Stats_summaryP <- Stats_summaryP %>% mutate(Percent_dNK3_of_Lymph = (dNK3/Total_lymph)*100)
Stats_summaryP <- Stats_summaryP %>% mutate(Percent_dNK1_of_NK = (dNK1/Total_NK)*100)
Stats_summaryP <- Stats_summaryP %>% mutate(Percent_dNK2_of_NK = (dNK2/Total_NK)*100)
Stats_summaryP <- Stats_summaryP %>% mutate(Percent_dNK3_of_NK = (dNK3/Total_NK)*100)

#Just copy and pasted the sheet into excel

#Identifying the differentially expressed genes in labour 

#Selecting the cell groups for decidua basalis

BdNK1_TNL <- Third %>% subset(subset = seurat_clusters == 5 & Condition == "TNL" & Tissue == "Basal Plate")
BdNK1_TNL <- rownames(BdNK1_TNL@meta.data)

BdNK1_TL <- Third %>% subset(subset = seurat_clusters == 5 & Condition == "TIL" & Tissue == "Basal Plate")
BdNK1_TL <- rownames(BdNK1_TL@meta.data)

BdNK2_TNL <- Third %>% subset(subset = seurat_clusters == 1 & Condition == "TNL" & Tissue == "Basal Plate")
BdNK2_TNL <- rownames(BdNK2_TNL@meta.data)

BdNK2_TL <- Third %>% subset(subset = seurat_clusters == 1 & Condition == "TIL" & Tissue == "Basal Plate")
BdNK2_TL <- rownames(BdNK2_TL@meta.data)

BdNK3_TNL <- Third %>% subset(subset = seurat_clusters == 2 & Condition == "TNL" & Tissue == "Basal Plate")
BdNK3_TNL <- rownames(BdNK3_TNL@meta.data)

BdNK3_TL <- Third %>% subset(subset = seurat_clusters == 2 & Condition == "TIL" & Tissue == "Basal Plate")
BdNK3_TL <- rownames(BdNK3_TL@meta.data)

#Selecting the cell groups for decidua parietalis 

dNK1_TNL <- Third %>% subset(subset = seurat_clusters == 5 & Condition == "TNL" & Tissue == "Chorioamniotic membranes")
dNK1_TNL <- rownames(dNK1_TNL@meta.data)

dNK1_TL <- Third %>% subset(subset = seurat_clusters == 5 & Condition == "TIL" & Tissue == "Chorioamniotic membranes")
dNK1_TL <- rownames(dNK1_TL@meta.data)

dNK2_TNL <- Third %>% subset(subset = seurat_clusters == 1 & Condition == "TNL" & Tissue == "Chorioamniotic membranes")
dNK2_TNL <- rownames(dNK2_TNL@meta.data)

dNK2_TL <- Third %>% subset(subset = seurat_clusters == 1 & Condition == "TIL" & Tissue == "Chorioamniotic membranes")
dNK2_TL <- rownames(dNK2_TL@meta.data)

dNK3_TNL <- Third %>% subset(subset = seurat_clusters == 2 & Condition == "TNL" & Tissue == "Chorioamniotic membranes")
dNK3_TNL <- rownames(dNK3_TNL@meta.data)

dNK3_TL <- Third %>% subset(subset = seurat_clusters == 2 & Condition == "TIL" & Tissue == "Chorioamniotic membranes")
dNK3_TL <- rownames(dNK3_TL@meta.data)

#Using find markers to compare TNL to TIL for each group, and find differentially expressed genes.
#Only select genes which have a p < 0.05

FindMarkers(Third, ident.1 = dNK1_TNL, ident.2 = dNK1_TL) %>% filter(p_val_adj <= 0.05) %>% write.csv(file = "dNK1_DP.csv")
FindMarkers(Third, ident.1 = dNK2_TNL, ident.2 = dNK2_TL) %>% filter(p_val_adj <= 0.05) %>% write.csv(file = "dNK2_DP.csv")
FindMarkers(Third, ident.1 = dNK3_TNL, ident.2 = dNK3_TL) %>% filter(p_val_adj <= 0.05) %>% write.csv(file = "dNK3_DP.csv")

FindMarkers(Third, ident.1 = BdNK1_TNL, ident.2 = BdNK1_TL) %>% filter(p_val_adj <= 0.05) %>% write.csv(file = "dNK1_DB.csv")
FindMarkers(Third, ident.1 = BdNK2_TNL, ident.2 = BdNK2_TL) %>% filter(p_val_adj <= 0.05) %>% write.csv(file = "dNK2_DB.csv")
FindMarkers(Third, ident.1 = BdNK3_TNL, ident.2 = BdNK3_TL) %>% filter(p_val_adj <= 0.05) %>% write.csv(file = "dNK3_DB.csv")




