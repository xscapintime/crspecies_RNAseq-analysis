# 2D pathway dot plot
# GSEA
# -------------------

rm(list = ls())

library(tidyverse)
library(ggplot2)
library(ggrepel)


## gsea
## ----
load("human_fgseaResTidy.Rdata")
load("mouse_fgseaResTidy.Rdata")

# overlap pathway
overlap_path <- intersect(h_fgseaResTidy$pathway, m_fgseaResTidy$pathway)

# select
human_pathway <- h_fgseaResTidy %>% filter(pathway %in% overlap_path) # padj < 0.05 &

mouse_pathway <- m_fgseaResTidy %>% filter(pathway %in% overlap_path) # padj < 0.05 &

# data to plot
# x as human, y as mouse
dat <- human_pathway[, c(1,3,6,7)]
dat <- inner_join(dat, mouse_pathway[, c(1,3,6,7)], by = c("pathway" = "pathway"))

dat$group <- ifelse(dat$NES.x * dat$NES.y > 0,
                    ifelse(abs(dat$NES.x) >= 1.5 & abs(dat$NES.y) >= 1.5 & dat$padj.x <= 0.05 & dat$padj.y <= 0.05, "corr", "other"),
                    ifelse(abs(dat$NES.x) >= 1.5 & abs(dat$NES.y) >= 1.5 & dat$padj.x <= 0.05 & dat$padj.y <= 0.05, "anti", "other"))

dat$sig <- ifelse(dat$padj.x <= 0.05 & dat$padj.y <= 0.05, "sig", "no")

# dotplot
theme_set(
  theme_classic() +
    theme(legend.position = "none"))


p <- ggplot(dat, #%>% filter(padj.x <= 0.05 | padj.y <= 0.05),
            aes(x = NES.x, y = NES.y))
p + geom_point(alpha = .6, aes(color = group, shape = sig,
                                size = dat$size.x/dat$size.y)) +

    scale_color_manual(values = c("#D43F3A", "#5CB85C", "#7C878E")) +

    # geom_text(data = dat %>%
    #           filter(group == "corr" &
    #           padj.x <= 0.01 & padj.y <= 0.01),
    #           mapping = aes(label = pathway), size = 3) +

    labs(title = "2D GSEA pathways",
          x = "Arrested / 8C NES", y = "Diapause / E4.5 NES")

ggsave(width = 7.6, height = 7.6, filename = "figs/2dgsea_nesall.png")



#

