# 2D pathway dot plot
# GSEA
# log2FC of pathway
# -------------------

rm(list = ls())

library(tidyverse)
library(ggplot2)
library(ggrepel)


## GSEA result
load("human_fgseaResTidy.Rdata")
load("mouse_fgseaResTidy.Rdata")


## pahtway logFC
human_pathway_fc <- read.table("human_pathway_fc.txt")
mouse_pathway_fc <- read.table("mouse_pathway_fc.txt")


## join
human_dat <- inner_join(h_fgseaResTidy, human_pathway_fc,
                        by = c("pathway" = "V1"))
colnames(human_dat)[9] <- "pathfc"

mouse_dat <- inner_join(m_fgseaResTidy, mouse_pathway_fc,
                        by = c("pathway" = "V1"))
colnames(mouse_dat)[9] <- "pathfc"

setdiff(human_dat$pathway, mouse_dat$pathway)


## overlap pathway
# overlap_path <- intersect(human_dat$pathway, mouse_dat$pathway)

# # select
# human_pathway <- human_dat %>% filter(pathway %in% overlap_path) # padj < 0.05 &

# mouse_pathway <- mouse_dat %>% filter(pathway %in% overlap_path) # padj < 0.05 &

# # check diff
# setdiff(human_pathway$pathway, mouse_pathway$pathway)



## data to plot
# x as human, y as mouse

dat_h <- human_dat %>%
        dplyr::select(pathway, pathfc, size, padj)

dat_m <- mouse_dat %>%
        dplyr::select(pathway, pathfc, size, padj)


dat <- inner_join(dat_h, dat_m, by = "pathway")


## clean up
dat <- dat %>% filter(size.x < 10) %>%
                filter(padj.x <= 0.1 | padj.y <= 0.1)




## plot data
# both fc and pval, no result
# dat$group <- ifelse(dat$pathfc.x * dat$pathfc.y > 0,
#                     ifelse(abs(dat$pathfc.x) >= 1 & abs(dat$pathfc.y) >= 1 & dat$padj.x <= 0.05 & dat$padj.y <= 0.05, "corr", "other"),
#                     ifelse(abs(dat$pathfc.x) >= 1 & abs(dat$pathfc.y) >= 1 & dat$padj.x <= 0.05 & dat$padj.y <= 0.05, "anti", "other"))


# only fc
dat$group <- ifelse(dat$pathfc.x * dat$pathfc.y > 0,
                    ifelse(abs(dat$pathfc.x) >= 1 & abs(dat$pathfc.y) >= 1, "corr", "other"),
                    ifelse(abs(dat$pathfc.x) >= 1 & abs(dat$pathfc.y) >= 1, "anti", "other"))


dat$sig <- ifelse(dat$padj.x <= 0.05 & dat$padj.y <= 0.05, "bothsig",
                    ifelse(dat$padj.x <= 0.05, "hsig",
                            ifelse(dat$padj.y <= 0.05, "msig", "no")))


# dotplot
theme_set(
  theme_classic() +
    theme(legend.position = "right"))


p <- ggplot(dat, #%>% filter(padj.x <= 0.05 | padj.y <= 0.05),
            aes(x = pathfc.x, y = pathfc.y))
p + geom_point(alpha = .6, aes(color = group, shape = sig)) +
                                #size = dat$size.x/dat$size.y)) + # size are the same

    scale_color_manual(values = c("#D43F3A", "#5CB85C", "#7C878E")) +

    # geom_text(data = dat %>%
    #           filter(group == "corr" &
    #           padj.x <= 0.01 & padj.y <= 0.01),
    #           mapping = aes(label = pathway), size = 3) +

    labs(title = "2D GSEA pathways",
          x = "Arrested / 8C log2FC", y = "Diapause / E4.5 log2FC")

ggsave(width = 7.6, height = 7.6, filename = "figs/2dgsea_pathfcsub.png")

