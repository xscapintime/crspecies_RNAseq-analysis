rm(list = ls())
#options(scipen = 999)
library(tidyverse)

library(DESeq2)

## read count (genes below cutoff removed)
# huamn
human_exp_mat <- read.table("human_select_rdcounts.tsv")

# mouse
mouse_exp_mat <- read.table("mouse_select_rdcounts.tsv")


## DEG
# create DESeq object
# design matrix
human_meta <- data.frame(stage = factor(c(rep("arr", each = 8),
                                        rep("8C", each = 20))))
row.names(human_meta) <- colnames(human_exp_mat)

mouse_meta <- data.frame(stage = factor(c(rep("dia", each = 2),
                                        rep("E4.5", each = 3))))
row.names(mouse_meta) <- colnames(mouse_exp_mat)


# Create DESeq2Dataset object
human_dds <- DESeqDataSetFromMatrix(countData = human_exp_mat, colData = human_meta, design = ~ stage)

mouse_dds <- DESeqDataSetFromMatrix(countData = mouse_exp_mat, colData = mouse_meta, design = ~ stage)


# filter
human_dds <- estimateSizeFactors(human_dds)
mouse_dds <- estimateSizeFactors(mouse_dds)


filterh <- rowSums(counts(human_dds, normalized = TRUE)) > 0
table(filter_h)
# filter_h
#  TRUE
# 15523
human_dds <- human_dds[filter_h, ]


filterm <- rowSums(counts(mouse_dds, normalized = TRUE)) > 0
table(filter_m)
# filter_m
#  TRUE
# 10783
mouse_dds <- mouse_dds[filter_m, ]


# only keep the homology/orthology? PAIRS
human_symbol <- read.table("human_idsyb.tsv")
mouse_symbol_tohuman <- read.table("mouse_idsyb_mapped2hugo.tsv")

overlap <- intersect(human_symbol$hgnc_symbol, mouse_symbol_tohuman$hsapiens_homolog_associated_gene_name)

overlap_h <- human_symbol %>% filter(hgnc_symbol %in% overlap)
human_dds <- human_dds[overlap_h$ensembl_gene_id, ]

overlap_m <- mouse_symbol_tohuman %>% filter(hsapiens_homolog_associated_gene_name %in% overlap)
mouse_dds <- mouse_dds[overlap_m$ensembl_gene_id, ]


# DEG
human_res <- results(DESeq(human_dds), contrast = c("stage", "arr", "8C"),
                    tidy = TRUE)
write.table(human_res, file = "human_deg.tsv", quote = F, sep = "\t")

mouse_res2 <- results(DESeq(mouse_dds), contrast = c("stage", "E4.5", "dia"),
                    tidy = TRUE, independentFiltering = F)
write.table(mouse_res, file = "mouse_deg.tsv", quote = F, sep = "\t")