# logFc for pathways
# GSEA
# ------------------

rm(list = ls())
library(tidyverse)

## GSEA results
load("human_fgseaResTidy.Rdata")
load("mouse_fgseaResTidy.Rdata")


## DE results
human_res <- read.table("human_deg.tsv")
mouse_res <- read.table("mouse_deg.tsv")

# need remove NA pval/qval genes
human_res <- human_res %>% na.omit()
mouse_res <- mouse_res %>% na.omit()

# NA pval in mouse data, make them be in the same length
human_res <- human_res %>% filter(row %in% mouse_res$row)


## pathway fold change
human_pathway_fc <- c()

system.time({
    
    for (i in seq_len(nrow(h_fgseaResTidy))) {

    pathway <- h_fgseaResTidy$pathway[i]

    de_inpath <- h_fgseaResTidy$leadingEdge[[i]]

    fc <- human_res %>%
            filter(row %in% de_inpath) %>%
            select(log2FoldChange) %>% summarise_all(mean)
    human_pathway_fc[i] <- fc
    names(human_pathway_fc)[i] <- pathway
    }
})

tmp <- unlist(human_pathway_fc)