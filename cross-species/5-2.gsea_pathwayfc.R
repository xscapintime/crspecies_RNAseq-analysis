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
    human_pathway_fc[i] <- fc %>% unlist()
    names(human_pathway_fc)[i] <- pathway

    }
})

write.table(human_pathway_fc, file = "human_pathway_fc.txt",
            quote = F, sep = "\t", row.names = T, col.names = F)


mouse_pathway_fc <- c()

system.time({
    
    for (i in seq_len(nrow(m_fgseaResTidy))) {

    pathway <- m_fgseaResTidy$pathway[i]

    de_inpath <- m_fgseaResTidy$leadingEdge[[i]]

    fc <- mouse_res %>%
            filter(row %in% de_inpath) %>%
            select(log2FoldChange) %>% summarise_all(mean)
    mouse_pathway_fc[i] <- fc %>% unlist()
    names(mouse_pathway_fc)[i] <- pathway

    }
})

write.table(mouse_pathway_fc, file = "mouse_pathway_fc.txt",
            quote = F, sep = "\t", row.names = T, col.names = F)




