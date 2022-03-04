#Used the fetch data tool to calculate the percentages of dNK1-3 from lymphocytes and from total dNK

#Loading libraries
library(Seurat)
library(patchwork)
library(dplyr)
library(tidyr)
library(openxlsx)

#setting working directory
setwd("*")

#Loading integrated data
all <- readRDS("Fourdata_paper.rds")

#Identified clusters in "analysing intergrated data" using a combo of first trimester "annotations" and cluster markers

#Cluster 3 = dNK1, 1 = dNK2, 5 = dNK3 and 6 = contains dNKp but also proliferating T cells. Therefore used GNLY expression to subset to just dNKp
DefaultAssay(all) <- "RNA"
stats <- FetchData(all, vars = c("Time", "Sample_Tag", "seurat_clusters", "GNLY"))

#remove as no NK cells in it so makes coding dNKP difficult
stats <- stats %>% filter(Sample_Tag != "NonPreg_SAM50")
stats <- stats %>% filter(Sample_Tag != "NonPreg_SAM46")

#Generating a thing to count
cell <- c(rep(1, times = 30489))
stats <- mutate(stats, cell)

#Adding all the individual columns
# The seurat clusters exlcuded for the Total lymph count are myeloid in origin

Stats_summary <- stats %>% 
  filter(seurat_clusters != 2, seurat_clusters != 4, seurat_clusters != 10, seurat_clusters != 11, seurat_clusters != 12) %>%
  group_by(Sample_Tag, Time) %>%
  arrange(Sample_Tag) %>%
  summarize(Total_lymph = sum(cell))

dNK1 <- stats %>% 
  group_by(Sample_Tag, seurat_clusters) %>%
  summarize(N = n()) %>%
  ungroup() %>%
  complete(Sample_Tag, seurat_clusters, fill = list(N=0)) %>%
  filter(seurat_clusters == 3) %>%
  arrange(Sample_Tag) %>%
  select(N)

Stats_summary$dNK1 <- dNK1

dNK2 <- stats %>% 
  group_by(Sample_Tag, seurat_clusters) %>%
  summarize(N = n()) %>%
  ungroup() %>%
  complete(Sample_Tag, seurat_clusters, fill = list(N=0)) %>%
  filter(seurat_clusters == 1) %>%
  arrange(Sample_Tag) %>%
  select(N)

Stats_summary$dNK2 <- dNK2

dNK3 <- stats %>% 
  group_by(Sample_Tag, seurat_clusters) %>%
  summarize(N = n()) %>%
  ungroup() %>%
  complete(Sample_Tag, seurat_clusters, fill = list(N=0)) %>%
  filter(seurat_clusters == 5) %>%
  arrange(Sample_Tag) %>%
  select(N)

Stats_summary$dNK3 <- dNK3

dNKp <- stats %>% 
  filter(GNLY >= 1) %>%
  group_by(Sample_Tag, seurat_clusters) %>%
  summarize(N = n()) %>%
  ungroup() %>%
  complete(Sample_Tag, seurat_clusters, fill = list(N=0)) %>%
  filter(seurat_clusters == 6) %>%
  arrange(Sample_Tag) %>%
  select(N)

Stats_summary$dNKp <- dNKp

Total_immune <- stats %>% 
  group_by(Sample_Tag, Time) %>%
  arrange(Sample_Tag) %>%
  summarize(N = sum(cell)) %>%
  ungroup() %>%
  complete(Sample_Tag, fill = list(N=0)) %>%
  arrange(Sample_Tag) %>%
  select(N)


Stats_summary$Total_Immune <- Total_immune

#Use mutate to calculate the percentages 

Stats_summary <- Stats_summary %>% mutate(Total_NK = sum(dNK1, dNK2, dNK3, dNKp))
Stats_summary <- Stats_summary %>% mutate(Percent_NK_of_Lymph = (Total_NK/Total_lymph)*100)
Stats_summary <- Stats_summary %>% mutate(Percent_dNK1_of_Lymph = (dNK1/Total_lymph)*100)
Stats_summary <- Stats_summary %>% mutate(Percent_dNK2_of_Lymph = (dNK2/Total_lymph)*100)
Stats_summary <- Stats_summary %>% mutate(Percent_dNK3_of_Lymph = (dNK3/Total_lymph)*100)
Stats_summary <- Stats_summary %>% mutate(Percent_dNK1_of_NK = (dNK1/Total_NK)*100)
Stats_summary <- Stats_summary %>% mutate(Percent_dNK2_of_NK = (dNK2/Total_NK)*100)
Stats_summary <- Stats_summary %>% mutate(Percent_dNK3_of_NK = (dNK3/Total_NK)*100)

#Just copy and pasted the sheet into excel





