rm(list = ls())
options(scipen = 999)
library(tidyverse)

library(DESeq2)

## data
# human
human_gene_mat <- read.table("human_readcounts.tsv")

# mouse
mouse_gene_mat <- read.table("mouse_readcounts.tsv")


## DESeq2 object
# convert to integer
#count_mx <- apply(count_mx, 2, as.integer)
#dimnames(count_mx) <- list(genes, samples)

# design matrix
human_meta <- data.frame(stage = factor(c(rep("arr", each = 8),
                                        rep("8C", each = 20))))
row.names(human_meta) <- colnames(human_gene_mat)

mouse_meta <- data.frame(stage = factor(c(rep("dia", each = 2),
                                        rep("E4.5", each = 3))))
row.names(mouse_meta) <- colnames(mouse_gene_mat)


# Check that sample names match in both files
all(colnames(human_gene_mat) %in% rownames(human_meta))
all(colnames(human_gene_mat) == rownames(human_meta))

all(colnames(mouse_gene_mat) %in% rownames(mouse_meta))
all(colnames(mouse_gene_mat) == rownames(mouse_meta))


# Create DESeq2Dataset object
human_dds <- DESeqDataSetFromMatrix(countData = human_gene_mat, colData = human_meta, design = ~ stage)

mouse_dds <- DESeqDataSetFromMatrix(countData = mouse_gene_mat, colData = mouse_meta, design = ~ stage)


## DESeq2 normalization
human_dds <- estimateSizeFactors(human_dds)
human_normalized_counts <- counts(human_dds, normalized = TRUE)

mouse_dds <- estimateSizeFactors(mouse_dds)
mouse_normalized_counts <- counts(mouse_dds, normalized = TRUE)



# maximum sample-group-averaged expression level ####
## @ramseyMethodCrossSpeciesVisualization2018
## intra-sample-group geometric mean (respectively) of the log2 counts
## inter-sample-group maximum

library(reshape2)
human_maxexp <- setNames(aggregate(value ~ gene,
                        data = aggregate(value ~ gene + stage,
                        data = merge(setNames(melt(log2(1 + as.matrix(human_normalized_counts))),
                        c("gene", "sample", "value")), human_meta, by.x = "sample", by.y = 0),
                        FUN = mean), FUN = max), c("gene", "max_exp"))

mouse_maxexp <- setNames(aggregate(value ~ gene,
                        data = aggregate(value ~ gene + stage,
                        data = merge(setNames(melt(log2(1 + as.matrix(mouse_normalized_counts))),
                        c("gene", "sample", "value")), mouse_meta, by.x = "sample", by.y = 0),
                        FUN = mean), FUN = max), c("gene", "max_exp"))


# minimum-expression-level cutoff using kernel density estimation ####
human_maxexp_density <- density(human_maxexp$max_exp)
plot(human_maxexp_density)
cutoff_human <- optimize(approxfun(human_maxexp_density$x, human_maxexp_density$y),
                        interval = c(1, quantile(human_maxexp_density$x)[2]))$minimum

mouse_maxexp_density <- density(mouse_maxexp$max_exp)
plot(mouse_maxexp_density)
cutoff_mouse <- optimize(approxfun(mouse_maxexp_density$x, mouse_maxexp_density$y),
                        interval = c(1, quantile(mouse_maxexp_density$x)[2]))$minimum

# remove genes below cutoff
human_maxexp_genes <- human_maxexp$gene[which(human_maxexp$max_exp >= cutoff_human)]
human_exp_mat <- human_gene_mat[human_maxexp_genes, ]
human_exp_norm_mat <- human_normalized_counts[human_maxexp_genes, ] #
write.table(human_exp_mat, file = "human_select_rdcounts.tsv", quote = F, sep = "\t")

mouse_maxexp_genes <- mouse_maxexp$gene[which(mouse_maxexp$max_exp >= cutoff_mouse)]
mouse_exp_mat <- mouse_gene_mat[mouse_maxexp_genes, ]
mouse_exp_norm_mat <- mouse_normalized_counts[mouse_maxexp_genes, ] #
write.table(mouse_exp_mat, file = "mouse_select_rdcounts.tsv", quote = F, sep = "\t")