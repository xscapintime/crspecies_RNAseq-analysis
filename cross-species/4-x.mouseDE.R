# to test mouse data: diapause vs E4.5
# ------------------------------------

rm(list = ls())
options(scipen = 999)
library(tidyverse)

library(DESeq2)

## read count, from te_counter, TEs excluded
mouse_gene_mat <- read.table("mouse_readcounts.tsv")

## remove 0 count genes
mouse_gene_mat <- mouse_gene_mat[which(rowSums(mouse_gene_mat) > 0),]


## DEG
# create DESeq object
# design matrix
mouse_meta <- data.frame(stage = factor(c(rep("dia", each = 2),
                                        rep("E4.5", each = 3))))
row.names(mouse_meta) <- colnames(mouse_gene_mat)


# Create DESeq2Dataset object
mouse_dds <- DESeqDataSetFromMatrix(countData = mouse_gene_mat, colData = mouse_meta, design = ~ stage)


# DEG
mouse_res <- results(DESeq(mouse_dds), contrast = c("stage", "E4.5", "dia"),
                    tidy = TRUE)
#write.table(mouse_res, file = "mouse_excdeg.tsv", quote = F, sep = "\t")



## looks like the difference between diapause and E4.5 is not that big?
## thus the padj ...