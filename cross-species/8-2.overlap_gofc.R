# 2D pathway dot plot
# GO
# -------------------

rm(list = ls())

library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(ggpmisc)



## pahtway logFC
human_bp_fc <- read.table("human_bp_fc.txt", row.names = 1)
human_mf_fc <- read.table("human_mf_fc.txt", row.names = 1)
human_cc_fc <- read.table("human_cc_fc.txt", row.names = 1)

mouse_bp_fc <- read.table("mouse_bp_fc.txt", row.names = 1)
mouse_mf_fc <- read.table("mouse_mf_fc.txt", row.names = 1)
mouse_cc_fc <- read.table("mouse_cc_fc.txt", row.names = 1)



## topGO result
load("human_GOallres.Rdata")
load("mouse_GOallres.Rdata")


# BP
human_bp <- cbind(human_allres[["BP"]], human_bp_fc)
mouse_bp <- cbind(mouse_allres[["BP"]], mouse_bp_fc)


# MF
human_mf <- cbind(human_allres[["MF"]], human_mf_fc)
mouse_mf <- cbind(mouse_allres[["MF"]], mouse_mf_fc)

# CC
human_cc <- cbind(human_allres[["CC"]], human_cc_fc)
mouse_cc <- cbind(mouse_allres[["CC"]], mouse_cc_fc)



## overlap pathway
# BP
dat <- inner_join(human_bp, mouse_bp, by = "GO.ID") %>%
            dplyr::select(GO.ID, Term.x, Annotated.x, classicFisher.x, V2.x,
                            Annotated.y, classicFisher.y, V2.y)



dat$group <- ifelse(dat$V2.x * dat$V2.y > 0,
                    ifelse(abs(dat$V2.x) >= 1 & abs(dat$V2.y) >= 1, "corr", "other"),
                    ifelse(abs(dat$V2.x) >= 1 & abs(dat$V2.y) >= 1, "anti", "other"))


dat$sig <- ifelse(dat$classicFisher.x <= 0.05 & dat$classicFisher.y <= 0.05, "bothsig",
                    ifelse(dat$classicFisher.x <= 0.05, "hsig",
                            ifelse(dat$classicFisher.y <= 0.05, "msig", "no")))


# dotplot
theme_set(
  theme_classic() +
    theme(legend.position = "right"))

formula <- y ~ x

p <- ggplot(dat, #%>% filter(padj.x <= 0.05 | padj.y <= 0.05),
            aes(x = V2.x, y = V2.y))
p + geom_point(alpha = .6, aes(color = group, shape = sig)) +
                                #size = dat$size.x/dat$size.y)) + # size are the same

    scale_color_manual(values = c("#D43F3A", "#5CB85C", "#7C878E")) +

    # geom_text(data = dat %>%
    #           filter(group == "corr" &
    #           padj.x <= 0.01 & padj.y <= 0.01),
    #           mapping = aes(label = pathway), size = 3) +

    labs(title = "2D GO BP pathways",
          x = "Arrested / 8C log2FC", y = "Diapause / E4.5 log2FC") +
    stat_cor(method = "pearson") +
    stat_poly_eq(aes(label = ..eq.label..),
        formula = formula, parse = TRUE, geom = "text", vjust = 2, hjust = 1)


# MF


# CC
human_cc <- read.table()