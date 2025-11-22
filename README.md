# reproductive_cycle
Code used for paper "Dynamic changes in uterine NK cell subset frequency and function over the menstrual cycle and pregnancy" 2022 - [Link to article](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2022.880438/full)

## Purpose
The purpose of this project was to combine scRNAseq data sets from different stages of the reproductive cycle and see how the immune system changes through the cycle. Four previously published datasets were used in this intergration. All analysis was done in R.

## Data sets
The non pregnant data and first trimester data are open access. The second trimester and third trimester data require discussion with the authors to use.

The non pregnant data is from "Mapping the temporal and spatial dynamics of the human endometrium in vivo and in vitro" 2021 -  [Link to article](https://www.nature.com/articles/s41588-021-00972-2).
The first trimester data is from "Single-cell reconstruction of the early maternal–fetal interface in humans" 2018 -  [Link to article](https://www.nature.com/articles/s41586-018-0698-6).
The second and third trimester data are from "Single cell transcriptional signatures of the human placenta in term and preterm parturition" 2019 - [Link to article](ttps://elifesciences.org/articles/52004)

## Subsetting 
There are four files explaining how each dataset was analysed and subset to just immune cells.
Subsetting non prenant to immune.R
Subsetting 1st trimester to immune.R
Subsetting 2nd trimester to immune.R
Subsetting 3rd trimester to immune.R

## Intergration
There is one script covering the intergration of these four immune datasets (Intergrate four datasets.R), one covering the analyse of the intergrated data (Analysing intergrated data.R).

## Visualisation
There are three scripts for visualising the intergrated data.
Calculating the frequency of uNK1-3.R
Generating violin plots for gene expression across the reproductive cycle.R
Calculating labouring data.R
