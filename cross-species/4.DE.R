rm(list = ls())
#options(scipen = 999)
library(tidyverse)

library(DESeq2)

## read count (genes below cutoff removed)
# huamn
human_exp_mat <- read.table("human_select_rdcounts.tsv")

# mouse
mouse_exp_mat <- read.table("mouse_select_rdcounts.tsv")



## only keep the homology/orthology? PAIRS
human_symbol <- read.table("human_idsyb.tsv")
mouse_symbol_tohuman <- read.table("mouse_idsyb_mapped2hugo.tsv")

index <- inner_join(human_symbol, mouse_symbol_tohuman, by = c("hgnc_symbol" = "hsapiens_homolog_associated_gene_name"))
write.table(index, file = "pairsidx.tsv", quote = F, sep = "\t")

# human count matrix
human_exp_mat$id <- row.names(human_exp_mat)
human_exp_mat <- inner_join(human_exp_mat, index[, 1:2], by = c("id" = "ensembl_gene_id.x"))

human_exp_mat <- human_exp_mat[, -29] %>% group_by(hgnc_symbol) %>% summarise_all(mean)

#human_exp_mat <- data.frame(human_exp_mat)
human_exp_mat$id <- index$ensembl_gene_id.x[match(human_exp_mat$hgnc_symbol, index$hgnc_symbol)]



# mouse count matrix
mouse_exp_mat$id <- row.names(mouse_exp_mat)
mouse_exp_mat <- inner_join(mouse_exp_mat, index[, 2:3], by = c("id" = "ensembl_gene_id.y"))

mouse_exp_mat <- mouse_exp_mat[, -6] %>% group_by(hgnc_symbol) %>% summarise_all(mean)

#mouse_exp_mat <- data.frame(mouse_exp_mat)
mouse_exp_mat$id <- index$ensembl_gene_id.y[match(mouse_exp_mat$hgnc_symbol, index$hgnc_symbol)]




## DEG
# create DESeq object
# design matrix
human_meta <- data.frame(stage = factor(c(rep("arr", each = 8),
                                        rep("8C", each = 20))))
row.names(human_meta) <- colnames(human_exp_mat)[2:29]

mouse_meta <- data.frame(stage = factor(c(rep("dia", each = 2),
                                        rep("E4.5", each = 3))))
row.names(mouse_meta) <- colnames(mouse_exp_mat)[2:6]


# Create DESeq2Dataset object
human_dat <- sapply(human_exp_mat[, 2:29], as.integer)

human_dds <- DESeqDataSetFromMatrix(countData = human_exp_mat[, 2:29], colData = human_meta, design = ~ stage)

mouse_dds <- DESeqDataSetFromMatrix(countData = mouse_exp_mat[, 2:6], colData = mouse_meta, design = ~ stage)


# filter
human_dds <- estimateSizeFactors(human_dds)
mouse_dds <- estimateSizeFactors(mouse_dds)


filter_h <- rowSums(counts(human_dds, normalized = TRUE)) > 0
table(filter_h)
# filter_h
#  TRUE
# 15523
human_dds <- human_dds[filter_h, ]


filter_m <- rowSums(counts(mouse_dds, normalized = TRUE)) > 0
table(filter_m)
# filter_m
#  TRUE
# 10783
mouse_dds <- mouse_dds[filter_m, ]



# DEG
human_res <- results(DESeq(human_dds), contrast = c("stage", "arr", "8C"),
                    tidy = TRUE)
write.table(human_res, file = "human_deg.tsv", quote = F, sep = "\t")

mouse_res2 <- results(DESeq(mouse_dds), contrast = c("stage", "E4.5", "dia"),
                    tidy = TRUE, independentFiltering = F)
write.table(mouse_res, file = "mouse_deg.tsv", quote = F, sep = "\t")