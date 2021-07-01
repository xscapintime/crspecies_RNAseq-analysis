# 2D pathway dot plot
# GO
# -------------------

rm(list = ls())

library(tidyverse)
library(ggplot2)
library(ggrepel)


## BP
human_bp <- read.table("human_predfres_bp.tsv", header = T, sep = "\t")
mouse_bp <- read.table("mouse_predfres_bp.tsv", header = T, sep = "\t")

# overlap pathway
overlap_path <- intersect(human_bp$GO.ID, mouse_bp$GO.ID)

# select
human_pathway <- h_fgseaResTidy %>% filter(pathway %in% overlap_path) # padj < 0.05 &

mouse_pathway <- m_fgseaResTidy %>% filter(pathway %in% overlap_path) # padj < 0.05 &

# MF


# CC
human_cc <- read.table()