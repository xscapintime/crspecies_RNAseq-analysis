rm(list = ls())
options(scipen = 999)
library(tidyverse)

library(DESeq2)
library(biomaRt)

## read count (genes below cutoff removed)
# huamn
human_exp_mat <- read.table("human_select_rdcounts.tsv")

# moue
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


# DEG
human_res <- results(DESeq(human_dds), tidy = TRUE)
write.table(human_res, file = "human_deg.tsv", quote = F, sep = "\t")

mouse_res <- results(DESeq(mouse_dds), tidy = TRUE)
write.table(mouse_res, file = "mouse_deg.tsv", quote = F, sep = "\t")


