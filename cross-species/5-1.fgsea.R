rm(list = ls())
#options(scipen = 999)
library(tidyverse)

## fgsea to do GSEA
#BiocManager::install("fgsea")
library(fgsea)


## DEG table
human_res <- read.table("human_deg.tsv")
mouse_res <- read.table("mouse_deg.tsv")

# need remove NA pval/qval genes
human_res <- human_res %>% na.omit()
mouse_res <- mouse_res %>% na.omit()

# NA pval in mouse data, make them be in the same length
human_res <- human_res %>% filter(row %in% mouse_res$row)


## ensembl id and symble
# human_symbol <- read.table("human_idsyb.tsv")
# mouse_symbol_tohuman <- read.table("mouse_idsyb_mapped2hugo.tsv")


# ## join deg table and symbol
# human_res_syb <- inner_join(human_res, human_symbol,
#                             by = c("row" = "ensembl_gene_id"))

# mouse_res_syb <- inner_join(mouse_res, mouse_symbol_tohuman,
#                             by = c("row" = "ensembl_gene_id"))


## ranked list (ranked by stat, Wald statistic)
human_stat <- human_res %>%
        dplyr::select(row, stat) %>%
        na.omit() %>%
        distinct() %>%
        group_by(row) %>%
        summarize(stat = mean(stat))
h_ranks <- deframe(human_stat)
barplot(sort(h_ranks, decreasing = T))


mouse_stat <- mouse_res %>%
        dplyr::select(row, stat) %>%
        na.omit() %>%
        distinct() %>%
        group_by(row) %>%
        summarize(stat = mean(stat))
m_ranks <- deframe(mouse_stat)
barplot(sort(m_ranks, decreasing = T))



### fgsea
### -----
## MSigDB 7.4 all gene sets
pathways.all <- gmtPathways("gmtdata/msigdb.v7.4.symbols.gmt")


## gsea
h_fgseaRes <- fgsea(pathways = pathways.all, h_ranks)
h_fgseaResTidy <- h_fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES))
save(h_fgseaResTidy, file = "human_fgseaResTidy.Rdata")

m_fgseaRes <- fgsea(pathways = pathways.all, m_ranks)
m_fgseaResTidy <- m_fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES))
save(m_fgseaResTidy, file = "mouse_fgseaResTidy.Rdata")