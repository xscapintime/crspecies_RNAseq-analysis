## base line expression as cutoff

rm(list = ls())
options(scipen = 999)
library(tidyverse)

library(DESeq2)

## data
# human
human_gene_mat_t2 <- read.table("human_readcounts_t2.tsv")
human_gene_mat_t3 <- read.table("human_readcounts_t3.tsv")

# mouse
mouse_gene_mat <- read.table("mouse_readcounts.tsv")


## DESeq2 object
# convert to integer
#count_mx <- apply(count_mx, 2, as.integer)
#dimnames(count_mx) <- list(genes, samples)

# design matrix
human_meta_t2 <- data.frame(stage = factor(c(rep("arr_t2", each = 5), rep("morula", each = 14))))
row.names(human_meta_t2) <- colnames(human_gene_mat_t2)

human_meta_t3 <- data.frame(stage = factor(c(rep("arr_t3", each = 8), rep("E4", each = 117))))
row.names(human_meta_t3) <- colnames(human_gene_mat_t3)


mouse_meta <- data.frame(stage = factor(c(rep("dia", each = 2),
                                        rep("E4.5", each = 3))))
row.names(mouse_meta) <- colnames(mouse_gene_mat)


# Check that sample names match in both files
all(colnames(human_gene_mat_t2) %in% rownames(human_meta_t2))
all(colnames(human_gene_mat_t2) == rownames(human_meta_t2))

all(colnames(human_gene_mat_t3) %in% rownames(human_meta_t3))
all(colnames(human_gene_mat_t3) == rownames(human_meta_t3))


all(colnames(mouse_gene_mat) %in% rownames(mouse_meta))
all(colnames(mouse_gene_mat) == rownames(mouse_meta))


# Create DESeq2Dataset object
human_dds_t2 <- DESeqDataSetFromMatrix(countData = human_gene_mat_t2, colData = human_meta_t2, design = ~ stage)
human_dds_t3 <- DESeqDataSetFromMatrix(countData = human_gene_mat_t3, colData = human_meta_t3, design = ~ stage)

mouse_dds <- DESeqDataSetFromMatrix(countData = mouse_gene_mat, colData = mouse_meta, design = ~ stage)


## DESeq2 normalization
human_dds_t2 <- estimateSizeFactors(human_dds_t2)
human_normalized_counts_t2 <- counts(human_dds_t2, normalized = TRUE)

human_dds_t3 <- estimateSizeFactors(human_dds_t3)
human_normalized_counts_t3 <- counts(human_dds_t3, normalized = TRUE)


mouse_dds <- estimateSizeFactors(mouse_dds)
mouse_normalized_counts <- counts(mouse_dds, normalized = TRUE)

# save normed count table
write.table(human_normalized_counts_t2, file = "human_normed_count_t2.tsv", quote = F, sep = "\t")
write.table(human_normalized_counts_t3, file = "human_normed_count_t3.tsv", quote = F, sep = "\t")

write.table(mouse_normalized_counts, file = "mouse_normed_count.tsv", quote = F, sep = "\t")


# maximum sample-group-averaged expression level ####
## @ramseyMethodCrossSpeciesVisualization2018
## intra-sample-group geometric mean (respectively) of the log2 counts
## inter-sample-group maximum

library(reshape2)
human_maxexp_t2 <- setNames(aggregate(value ~ gene,
                        data = aggregate(value ~ gene + stage,
                        data = merge(setNames(melt(log2(1 + as.matrix(human_normalized_counts_t2))),
                        c("gene", "sample", "value")), human_meta_t2, by.x = "sample", by.y = 0),
                        FUN = mean), FUN = max), c("gene", "max_exp"))

human_maxexp_t3 <- setNames(aggregate(value ~ gene,
                        data = aggregate(value ~ gene + stage,
                        data = merge(setNames(melt(log2(1 + as.matrix(human_normalized_counts_t3))),
                        c("gene", "sample", "value")), human_meta_t3, by.x = "sample", by.y = 0),
                        FUN = mean), FUN = max), c("gene", "max_exp"))


mouse_maxexp <- setNames(aggregate(value ~ gene,
                        data = aggregate(value ~ gene + stage,
                        data = merge(setNames(melt(log2(1 + as.matrix(mouse_normalized_counts))),
                        c("gene", "sample", "value")), mouse_meta, by.x = "sample", by.y = 0),
                        FUN = mean), FUN = max), c("gene", "max_exp"))


# minimum-expression-level cutoff using kernel density estimation ####
human_maxexp_density_t2 <- density(human_maxexp_t2$max_exp)
plot(human_maxexp_density_t2)
cutoff_human_t2 <- optimize(approxfun(human_maxexp_density_t2$x, human_maxexp_density_t2$y),
                        interval = c(1, quantile(human_maxexp_density_t2$x)[2]))$minimum

human_maxexp_density_t3 <- density(human_maxexp_t3$max_exp)
plot(human_maxexp_density_t3)
cutoff_human_t3 <- optimize(approxfun(human_maxexp_density_t3$x, human_maxexp_density_t3$y),
                        interval = c(1, quantile(human_maxexp_density_t3$x)[2]))$minimum


mouse_maxexp_density <- density(mouse_maxexp$max_exp)
plot(mouse_maxexp_density)
cutoff_mouse <- optimize(approxfun(mouse_maxexp_density$x, mouse_maxexp_density$y),
                        interval = c(1, quantile(mouse_maxexp_density$x)[2]))$minimum

# density plot together
library(ggplot2)

dat <- rbind(human_maxexp_t2, mouse_maxexp) # combine type 2 and type 3 with mouse data

dat$species <- factor(c(rep("human", each = nrow(human_maxexp_t2)), # type 2 or type 3
                        rep("mouse", each = nrow(mouse_maxexp))))

theme_set(
    theme_bw() +
    theme(legend.position = "right")
)

p <- ggplot(dat, aes(x = max_exp))
p + geom_density(aes(group = species, fill = species), size = .6, alpha = 0.6) +
    scale_color_manual(values = c("#1a7bca", "#EFC000FF")) +
    scale_fill_manual(values = c("#1a7bca", "#EFC000FF")) +
    geom_vline(aes(xintercept = cutoff_human_t2),  ## human
                linetype = "dashed", size = 0.6, color = "#1a7bca") +
    geom_text(aes(x = cutoff_human_t2, y = 0.15,
            label = paste0("Human cutoff: ", round(cutoff_human_t2, 2)))) +
    geom_vline(aes(xintercept = cutoff_mouse),  ## human
                linetype = "dashed", size = 0.6, color = "#EFC000FF") +
    geom_text(aes(x = cutoff_mouse, y = 0.1,
            label = paste0("Mouse cutoff: ", round(cutoff_mouse, 2)))) +
    xlab("Max normalized log2 expression")

ggsave(filename = "../figs/type2_Max_normalized_log2expression_dens.png")


### type 3
dat <- rbind(human_maxexp_t3, mouse_maxexp) # combine type 2 and type 3 with mouse data

dat$species <- factor(c(rep("human", each = nrow(human_maxexp_t3)), # type 2 or type 3
                        rep("mouse", each = nrow(mouse_maxexp))))

theme_set(
    theme_bw() +
    theme(legend.position = "right")
)

p <- ggplot(dat, aes(x = max_exp))
p + geom_density(aes(group = species, fill = species), size = .6, alpha = 0.6) +
    scale_color_manual(values = c("#1a7bca", "#EFC000FF")) +
    scale_fill_manual(values = c("#1a7bca", "#EFC000FF")) +
    geom_vline(aes(xintercept = cutoff_human_t3),  ## human
                linetype = "dashed", size = 0.6, color = "#1a7bca") +
    geom_text(aes(x = cutoff_human_t3, y = 0.15,
            label = paste0("Human cutoff: ", round(cutoff_human_t3, 2)))) +
    geom_vline(aes(xintercept = cutoff_mouse),  ## human
                linetype = "dashed", size = 0.6, color = "#EFC000FF") +
    geom_text(aes(x = cutoff_mouse, y = 0.1,
            label = paste0("Mouse cutoff: ", round(cutoff_mouse, 2)))) +
    xlab("Max normalized log2 expression")

ggsave(filename = "../figs/type3_Max_normalized_log2expression_dens.png")




# remove genes below cutoff
human_maxexp_genes_t2 <- human_maxexp_t2$gene[which(human_maxexp_t2$max_exp >= cutoff_human_t2)]
human_exp_mat_t2 <- human_gene_mat_t2[human_maxexp_genes_t2, ]
human_exp_norm_mat_t2 <- human_normalized_counts_t2[human_maxexp_genes_t2, ]
write.table(human_exp_mat_t2, file = "human_select_rdcounts_t2.tsv", quote = F, sep = "\t")

human_maxexp_genes_t3 <- human_maxexp_t3$gene[which(human_maxexp_t3$max_exp >= cutoff_human_t3)]
human_exp_mat_t3 <- human_gene_mat_t3[human_maxexp_genes_t3, ]
human_exp_norm_mat_t3 <- human_normalized_counts_t3[human_maxexp_genes_t3, ]
write.table(human_exp_mat_t3, file = "human_select_rdcounts_t3.tsv", quote = F, sep = "\t")


## symbolic link to file generated last time
# mouse_maxexp_genes <- mouse_maxexp$gene[which(mouse_maxexp$max_exp >= cutoff_mouse)]
# mouse_exp_mat <- mouse_gene_mat[mouse_maxexp_genes, ]
# mouse_exp_norm_mat <- mouse_normalized_counts[mouse_maxexp_genes, ] #
# write.table(mouse_exp_mat, file = "mouse_select_rdcounts.tsv", quote = F, sep = "\t")