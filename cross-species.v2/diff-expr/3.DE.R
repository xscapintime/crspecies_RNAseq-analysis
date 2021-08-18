## diffrential expresison analysis
## human ones


rm(list = ls())
#options(scipen = 999)
library(tidyverse)

library(DESeq2)

## read count (genes below cutoff removed)
# huamn
human_exp_mat_t2 <- read.table("human_select_rdcounts_t2.tsv")
human_exp_mat_t3 <- read.table("human_select_rdcounts_t3.tsv")

# mouse
mouse_exp_mat <- read.table("mouse_select_rdcounts.tsv")


## only keep the homology/orthology? PAIRS
mouse_symbol_tohuman <- read.table("mouse_idsyb_mapped2hugo.tsv")

# for type 2
human_symbol <- row.names(human_exp_mat_t2)
#index <- intersect(human_symbol, mouse_symbol_tohuman$hsapiens_homolog_associated_gene_name)
index_t2 <- inner_join(as.data.frame(human_symbol), mouse_symbol_tohuman, by = c("human_symbol" = "hsapiens_homolog_associated_gene_name"))
write.table(index_t2, file = "pairsidx_t2.tsv", quote = F, sep = "\t")

# for type 3
human_symbol <- row.names(human_exp_mat_t3)

#index <- intersect(human_symbol, mouse_symbol_tohuman$hsapiens_homolog_associated_gene_name)
index_t3 <- inner_join(as.data.frame(human_symbol), mouse_symbol_tohuman, by = c("human_symbol" = "hsapiens_homolog_associated_gene_name"))
write.table(index_t3, file = "pairsidx_t3.tsv", quote = F, sep = "\t")


## human count matrix
# typ2 2 data
human_exp_mat_t2 <- human_exp_mat_t2[unique(index_t2$human_symbol), ]
human_exp_mat_t3 <- human_exp_mat_t3[unique(index_t3$human_symbol), ]


## mouse count matrix
mouse_exp_mat$id <- row.names(mouse_exp_mat)

# for type 2 data
mouse_exp_mat_t2 <- inner_join(mouse_exp_mat, index_t2, by = c("id" = "ensembl_gene_id"))

mouse_exp_mat_t2 <- mouse_exp_mat_t2[, -6] %>% group_by(human_symbol) %>% summarise_all(mean)

#mouse_exp_mat <- data.frame(mouse_exp_mat)
mouse_exp_mat_t2$id <- index_t2$ensembl_gene_id[match(mouse_exp_mat_t2$human_symbol, index_t2$human_symbol)]


# for type 3 data
mouse_exp_mat_t3 <- inner_join(mouse_exp_mat, index_t3, by = c("id" = "ensembl_gene_id"))

mouse_exp_mat_t3 <- mouse_exp_mat_t3[, -6] %>% group_by(human_symbol) %>% summarise_all(mean)

#mouse_exp_mat <- data.frame(mouse_exp_mat)
mouse_exp_mat_t3$id <- index_t3$ensembl_gene_id[match(mouse_exp_mat_t3$human_symbol, index_t3$human_symbol)]



## DEG
# create DESeq object
# design matrix
human_meta_t2 <- data.frame(stage = factor(c(rep("arr_t2", each = 5), rep("morula", each = 14))))
row.names(human_meta_t2) <- colnames(human_exp_mat_t2)

human_meta_t3 <- data.frame(stage = factor(c(rep("arr_t3", each = 8), rep("E4", each = 117))))
row.names(human_meta_t3) <- colnames(human_exp_mat_t3)


mouse_meta <- data.frame(stage = factor(c(rep("dia", each = 2), rep("E4.5", each = 3))))
row.names(mouse_meta) <- colnames(mouse_exp_mat)[1:5]


# Create DESeq2Dataset object
human_dat_t2 <- sapply(human_exp_mat_t2, as.integer)
row.names(human_dat_t2) <- row.names(human_exp_mat_t2)

human_dat_t3 <- sapply(human_exp_mat_t3, as.integer)
row.names(human_dat_t3) <- row.names(human_exp_mat_t3)


mouse_dat_t2 <- sapply(mouse_exp_mat_t2[, 2:6], as.integer)
row.names(mouse_dat_t2) <- mouse_exp_mat_t2$human_symbol

mouse_dat_t3 <- sapply(mouse_exp_mat_t3[, 2:6], as.integer)
row.names(mouse_dat_t3) <- mouse_exp_mat_t3$human_symbol


human_dds_t2 <- DESeqDataSetFromMatrix(countData = human_dat_t2, colData = human_meta_t2, design = ~ stage)
human_dds_t3 <- DESeqDataSetFromMatrix(countData = human_dat_t3, colData = human_meta_t3, design = ~ stage)


mouse_dds_t2 <- DESeqDataSetFromMatrix(countData = mouse_dat_t2, colData = mouse_meta, design = ~ stage)
mouse_dds_t3 <- DESeqDataSetFromMatrix(countData = mouse_dat_t3, colData = mouse_meta, design = ~ stage)


# filter
human_dds_t2 <- estimateSizeFactors(human_dds_t2)
human_dds_t3 <- estimateSizeFactors(human_dds_t3)

mouse_dds_t2 <- estimateSizeFactors(mouse_dds_t2)
mouse_dds_t3 <- estimateSizeFactors(mouse_dds_t3)


filter_h <- rowSums(counts(mouse_dds_t2, normalized = TRUE)) > 0
table(filter_h)
# human_dds <- human_dds[filter_h, ]


filter_m <- rowSums(counts(mouse_dds_t3, normalized = TRUE)) > 0
table(filter_m)
# mouse_dds <- mouse_dds[filter_m, ]


### type 2 data 8103 genes
### type 3 data 8139 genes


# DEG
human_res_t2 <- results(DESeq(human_dds_t2), contrast = c("stage", "arr_t2", "morula"),
                    tidy = TRUE, independentFiltering = F)
write.table(human_res_t2, file = "human_deg_t2.tsv", quote = F, sep = "\t")

human_res_t3 <- results(DESeq(human_dds_t3), contrast = c("stage", "arr_t3", "E4"),
                    tidy = TRUE, independentFiltering = F)
write.table(human_res_t3, file = "human_deg_t3.tsv", quote = F, sep = "\t")



mouse_res_t2 <- results(DESeq(mouse_dds_t2), contrast = c("stage", "dia", "E4.5"),
                    tidy = TRUE, independentFiltering = F)
write.table(mouse_res_t2, file = "mouse_deg_t2.tsv", quote = F, sep = "\t")

mouse_res_t3 <- results(DESeq(mouse_dds_t3), contrast = c("stage", "dia", "E4.5"),
                    tidy = TRUE, independentFiltering = F)
write.table(mouse_res_t3, file = "mouse_deg_t3.tsv", quote = F, sep = "\t")