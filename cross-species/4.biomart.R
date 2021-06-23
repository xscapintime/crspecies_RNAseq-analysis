rm(list = ls())
options(scipen = 999)
library(tidyverse)

library(biomaRt)

## read count (genes below cutoff removed)
# huamn
human_exp_mat <- read.table("human_select_rdcounts.tsv")

# moue
mouse_exp_mat <- read.table("mouse_select_rdcounts.tsv")


## set ensembl datasets
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

## list attributes
listAttributes(human) %>% head


## ensembl id to hugo symbol
human_id <- row.names(human_exp_mat)
human_symbol <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                      filters = "ensembl_gene_id",
                      values = human_id, mart = human)
human_symbol  <- human_symbol %>% na_if("") %>% na.omit() %>% distinct()
write.table(human_symbol, file = "human_idsyb.tsv", quote = F, sep = "\t")

## ortholog pairs
mouse_id <- row.names(mouse_exp_mat)
mouse_symbol_tohuman <- getBM(attributes = c("ensembl_gene_id", "hsapiens_homolog_associated_gene_name"),
                              filters = "ensembl_gene_id",
                              values = mouse_id, mart = mouse)
mouse_symbol_tohuman  <- mouse_symbol_tohuman %>% na_if("") %>% na.omit() %>% distinct()
write.table(mouse_symbol_tohuman, file = "mouse_idsyb_mapped2hugo.tsv", quote = F, sep = "\t")



## no use
mouse_symbol_human <- getLDS(attributes = c("hgnc_symbol"),
                             filters = "hgnc_symbol", values = human_symbol$hgnc_symbol, mart = human,
                             attributesL = c("mgi_symbol", "ensembl_gene_id"), martL = mouse)