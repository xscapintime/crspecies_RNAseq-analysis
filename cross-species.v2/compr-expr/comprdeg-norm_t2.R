# compare read count of diapuse DE genes
# for type 2 arrested embryo data
# type 2 vs morula
# normalized read count, seleceted 8103 genes
# ------------------------------------------


rm(list = ls())

library(tidyverse)


## diapause DE gene set
load("../compr2d/diapuse_set_t2.Rdata")

# all up genes
up <- union(diapuse_set_t2[["Diapause_up_t2"]], diapuse_set_t2[["Diapause_up_Duy2021"]])
dn <- union(diapuse_set_t2[["Diapause_dn_t2"]], diapuse_set_t2[["Diapause_dn_Duy2021"]])


## normed read count
human_exp_normed <- read.table("../diff-expr/human_normed_count_t2.tsv")
mouse_exp_normed <- read.table("../diff-expr/mouse_normed_count.tsv")


## mouse ensemble id to hugo
index <- read.table("../diff-expr/pairsidx_t2.tsv")


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


## up gene
up <- up[up %in% row.names(human_exp_normed)]

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
        ylab("Normalized read count") + xlab("")
    ggsave(width = 7.6, height = 7.6, filename = paste0(names(plot_dat[i]), "_normed_up_t2.png"), path = "../figs/t2_boxplot/normed-deg_boxplot")

}


## down gene
dn <- dn[dn %in% row.names(human_exp_normed)]

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
        ylab("Normalized read count") + xlab("")

    ggsave(width = 7.6, height = 7.6, filename = paste0(names(plot_dat[i]), "_normed_dn_t2.png"), path = "../figs/t2_boxplot/normed-deg_boxplot")
}
