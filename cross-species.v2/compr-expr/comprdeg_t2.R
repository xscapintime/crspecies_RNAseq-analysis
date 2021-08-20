# compare read count of diapuse DE genes
# for type 2 arrested embryo data
# type 2 vs morula
# original read count, seleceted 8103 genes
# ------------------------------------------


rm(list = ls())

library(tidyverse)


## diapause DE gene set
load("../compr2d/diapuse_set_t2.Rdata")

# all up genes
up <- union(diapuse_set_t2[["Diapause_up_t2"]], diapuse_set_t2[["Diapause_up_Duy2021"]])
dn <- union(diapuse_set_t2[["Diapause_dn_t2"]], diapuse_set_t2[["Diapause_dn_Duy2021"]])


## read count table
human_exp_mat <- read.table("../diff-expr/human_select_rdcounts_t2.tsv")
mouse_exp_mat <- read.table("../diff-expr/mouse_select_rdcounts.tsv")

## mouse ensemble id to hugo
index <- read.table("../diff-expr/pairsidx_t2.tsv")


human_exp_mat <- human_exp_mat[unique(index$human_symbol), ]

mouse_exp_mat$id <- row.names(mouse_exp_mat)
mouse_exp_mat <- inner_join(mouse_exp_mat, index, by = c("id" = "ensembl_gene_id"))
mouse_exp_mat <- mouse_exp_mat[, -6] %>% group_by(human_symbol) %>% summarise_all(mean)


## every DE gene in diapause sets
## boxplot

library(reshape2)

human <- melt(as.matrix(human_exp_mat))
human$stage <- ifelse(grepl("Arrested", human$Var2), "arr_t2", "morula")
colnames(human)[1:2] <- c("hugo", "phe")

mouse <- melt(mouse_exp_mat)
mouse$stage <- ifelse(grepl("DIA", mouse$variable), "dia", "E4.5")
colnames(mouse)[1:2] <- c("hugo", "phe")

setdiff(mouse$hugo %>% unique, human$hugo %>% unique)


## buil a plot data list obj
## up gene
up <- up[up %in% row.names(human_exp_mat)]

plot_dat <- list()

for(g in up) {

    plot_dat[[g]] <- rbind(human %>% filter(hugo == g), mouse %>% filter(hugo == g))
    plot_dat[[g]]$stage <- factor(plot_dat[[g]]$stage, levels = c("arr_t2", "morula", "dia", "E4.5"))
    levels(plot_dat[[g]]$stage) <- c("Type Ⅱ Arrested", "Morula", "Diapause", "E4.5")

}


## boxplot

library(ggplot2)
library(hrbrthemes)
library(viridis)

theme_set(
  theme_ipsum() +
  theme(
      legend.position = "none",
      plot.title = element_text(size = 11)
    ))

for(i in 1:length(plot_dat)) {
  p <- ggplot(plot_dat[[i]], aes(x=stage, y=value, fill=stage))
  p +  geom_boxplot() +
    scale_fill_viridis(discrete = TRUE, alpha = 0.6) +
    geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
    labs(title = names(plot_dat)[i], caption = "Diapuse up DEG") +
    ylab("Read count") + xlab("")

    ggsave(width = 7.6, height = 7.6, filename = paste0(names(plot_dat[i]), "_up_t2.png"), path = "../figs/t2_boxplot/deg_boxplot")
}



## down gene
dn <- dn[dn %in% row.names(human_exp_mat)]

plot_dat <- list()

for(g in dn) {

    plot_dat[[g]] <- rbind(human %>% filter(hugo == g), mouse %>% filter(hugo == g))
    plot_dat[[g]]$stage <- factor(plot_dat[[g]]$stage, levels = c("arr_t2", "morula", "dia", "E4.5"))
    levels(plot_dat[[g]]$stage) <- c("Type Ⅱ Arrested", "Morula", "Diapause", "E4.5")

}


## boxplot

for(i in 1:length(plot_dat)) {
  p <- ggplot(plot_dat[[i]], aes(x=stage, y=value, fill=stage))
  p +  geom_boxplot() +
    scale_fill_viridis(discrete = TRUE, alpha = 0.6) +
    geom_jitter(color = "black", size = 0.4, alpha = 0.9) +
    labs(title = names(plot_dat)[i], caption = "Diapuse down DEG") +
    ylab("Read count") + xlab("")

    ggsave(width = 7.6, height = 7.6, filename = paste0(names(plot_dat[i]), "_dn_t2.png"), path = "../figs/t2_boxplot/deg_boxplot")
}
