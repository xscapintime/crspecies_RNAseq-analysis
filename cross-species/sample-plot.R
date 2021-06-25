# sample distance heatmap and PCA
# -------------------------------

rm(list = ls())
options(scipen = 999)
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


# DESeq
human_dds2 <- DESeq(human_dds)

mouse_dds2 <- DESeq(mouse_dds)


## sample distance heatmap
human_vsd <- vst(human_dds2, blind = FALSE)

mouse_vsd <- vst(mouse_dds2, blind = FALSE)



# heatmap
# -------
# human
sampleDists <- dist(t(assay(human_vsd)))

library("RColorBrewer")
library(pheatmap)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(human_exp_mat)
colnames(sampleDistMatrix) <- NULL

pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         cellwidth = 12,
         cellheight = 12,
         filename = "figs/human_sampledist.png")

# mouse
sampleDists <- dist(t(assay(mouse_vsd)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(mouse_exp_mat)
colnames(sampleDistMatrix) <- NULL

pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         cellwidth = NA,
         cellheight = NA,
         filename = "figs/mouse_sampledist.png")



# PCA
# ---

# human
pcaData <- plotPCA(human_vsd, intgroup = "stage", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color = stage)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  ggsave("figs/human_pca.png")


# mouse
pcaData <- plotPCA(mouse_vsd, intgroup = "stage", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color = stage)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  ggsave("figs/mouse_pca.png")

