rm(list = ls())
#options(scipen = 999)
library(tidyverse)
library(reshape2)

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

human_exp_mat$id <- row.names(human_exp_mat)
human_exp_mat <- inner_join(human_exp_mat, index, by = c("id" = "ensembl_gene_id.x"))

tolong <- melt(human_exp_mat)

#tapply(human_exp_mat[, 1:28], INDEX = human_exp_mat$hgnc_symbol, FUN = mean)

tolong <- tolong %>%
        dplyr::select(id, hgnc_symbol, variable, value) %>%
        na.omit() %>%
        distinct() %>%
        group_by(hgnc_symbol, variable, id) %>%
        summarize(value = mean(value))

test <- dcast(tolong, id ~ variable, value.var = "value")

human_exp_mat <- human_exp_mat[index$ensembl_gene_id.x, ]

overlap_m <- mouse_symbol_tohuman %>% filter(hsapiens_homolog_associated_gene_name %in% overlap)
mouse_exp_mat <- mouse_exp_mat[overlap_m$ensembl_gene_id, ]



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