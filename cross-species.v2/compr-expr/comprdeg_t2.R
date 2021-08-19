# compare read count of diapuse DE genes
# for type 2 arrested embryo data
# type 2 vs morula

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



# buil a plot data list obj
# up gene
plot_dat <- list()

for(g in up) {
    plot_dat[[g]] <- rbind(human %>% filter(hugo == g), mouse %>% filter(hugo == g))
    plot_dat[[g]]$stage <- factor(plot_dat[[g]]$stage, levels = c("arr_t2", "morula", "dia", "E4.5"))

}



library(ggplot2)
library(hrbrthemes)
library(viridis)

ggplot(plot_dat[[40]], aes(x=stage, y=value, fill=stage)) +
    geom_boxplot() +
    scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    geom_jitter(color="black", size=0.4, alpha=0.9) +
    theme_ipsum() +
    theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) +
    ggtitle("A boxplot with jitter") +
    xlab("")




### if normalized?
human_exp_normed <- read.table("../diff-expr/human_normed_count_t2.tsv")
mouse_exp_normed <- read.table("../diff-expr/mouse_normed_count.tsv")

human_exp_normed <- human_exp_normed[unique(index$human_symbol), ]

mouse_exp_normed$id <- row.names(mouse_exp_normed)
mouse_exp_normed <- inner_join(mouse_exp_normed, index, by = c("id" = "ensembl_gene_id"))
mouse_exp_normed <- mouse_exp_normed[, -6] %>% group_by(human_symbol) %>% summarise_all(mean)


## plot data
human <- melt(as.matrix(human_exp_normed))
human$stage <- ifelse(grepl("Arrested", human$Var2), "arr_t2", "morula")
colnames(human)[1:2] <- c("hugo", "phe")

mouse <- melt(mouse_exp_normed)
mouse$stage <- ifelse(grepl("DIA", mouse$variable), "dia", "E4.5")
colnames(mouse)[1:2] <- c("hugo", "phe")


plot_dat <- list()

for(g in up) {
    plot_dat[[g]] <- rbind(human %>% filter(hugo == g), mouse %>% filter(hugo == g))
    plot_dat[[g]]$stage <- factor(plot_dat[[g]]$stage, levels = c("arr_t2", "morula", "dia", "E4.5"))

}

ggplot(plot_dat[[40]], aes(x=stage, y=value, fill=stage)) +
    geom_boxplot() +
    scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    geom_jitter(color="black", size=0.4, alpha=0.9) +
    theme_ipsum() +
    theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) +
    ggtitle("A boxplot with jitter") +
    xlab("")