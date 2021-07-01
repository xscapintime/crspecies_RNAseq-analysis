# logFc for pathways
# ORA GO
# ------------------

rm(list = ls())
library(tidyverse)
library(topGO)


## Godata to retrive genes in GO term
load("human_GOdata.Rdata")
load("mouse_GOdata.Rdata")

load("human_GOallres.Rdata")
load("mouse_GOallres.Rdata")


## DE results
human_res <- read.table("human_deg.tsv")
mouse_res <- read.table("mouse_deg.tsv")

# need remove NA pval/qval genes
human_res <- human_res %>% na.omit()
mouse_res <- mouse_res %>% na.omit()

# NA pval in mouse data, make them be in the same length
human_res <- human_res %>% filter(row %in% mouse_res$row)



## BP
human_bp <- human_allres[["BP"]]
mouse_bp <- mouse_allres[["BP"]]



selcGenes <- genesInTerm(myGOdata, whichGO=selcTerm)


