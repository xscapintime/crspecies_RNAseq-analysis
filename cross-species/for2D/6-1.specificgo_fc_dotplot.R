# all gene logFC dotplot
# x: mouse logFC, y: human logFC
# ------------------------------


rm(list = ls())

library(tidyverse)


## DESeq2 results
human_res <- read.table("../human_deg.tsv")
mouse_res <- read.table("../mouse_deg.tsv")

# need remove NA pval/qval genes
human_res <- human_res %>% na.omit()
mouse_res <- mouse_res %>% na.omit()

# NA pval in mouse data, make them be in the same length
human_res <- human_res %>% filter(row %in% mouse_res$row)



## genes mapped to gene sets
load("gene2sets_10cut.Rdata")


## pathways to look at, from Scognamiglio et al., 2016

# downs
dnidx <- c(names(map_fin) %>% grep(pattern = "GOBP_RNA_PROCESSING$"),
            names(map_fin) %>% grep(pattern = "GOBP_RNA_SPLICING$"),
            names(map_fin) %>% grep(pattern = "GOCC_SPLICEOSOMAL_COMPLEX$"),
            names(map_fin) %>% grep(pattern = "GOBP_RIBOSOME_BIOGENESIS$"),
            names(map_fin) %>% grep(pattern = "GOBP_DNA_REPLICATION$"))

# GO:0006412 not in MSigdb, translation


# ups
upidx <- c(names(map_fin) %>% grep(pattern = "GOMF_INTEGRIN_BINDING$"),
            names(map_fin) %>% grep(pattern = "GOMF_G_PROTEIN_COUPLED_RECEPTOR_ACTIVITY$"))

# GO:0038023 not in MSigdb, signaling receptor activity



### plot data
any(human_res$row != mouse_res$row)

dat <- cbind(human_res$log2FoldChange, mouse_res$log2FoldChange)
dimnames(dat) <- list(human_res$row, c("human_fc", "mouse_fc"))
dat <- data.frame(dat)



### figure setting
library(ggplot2)
library(ggthemes)
library(plyr)

theme_set(theme_few() + theme(legend.position = "right"))


####################################################
## for down pathways in Scognamiglio et al., 2016 ##
####################################################

plotdat <- list()
for (i in 1:length(dnidx)) {
    id <- dnidx[i]
    plotdat[[i]] <- cbind(dat, map_fin[[id]])
    colnames(plotdat[[i]])[3] <- "inpath"
    plotdat[[i]]$inpath <- factor(plotdat[[i]]$inpath)
    names(plotdat)[i] <- names(map_fin[id])
}


for (i in 1:length(plotdat)) {

    p <- ggplot()
    p + geom_point(data = plotdat[[i]] %>% filter(inpath == 0),
                    color = "#aab0b4", alpha = .6,
                    aes(x = mouse_fc, y = human_fc)) +
        geom_point(data = plotdat[[i]] %>% filter(inpath == 1),
                    color = "#186eb4", alpha = .6,
                    aes(x = mouse_fc, y = human_fc)) +

        #scale_color_manual(values = c("#aab0b4", "#186eb4")) +

        geom_vline(xintercept = 0, linetype = "dashed", size = 0.6) +
        geom_hline(yintercept = 0, linetype = "dashed", size = 0.6) +

        labs(title = names(plotdat[i]),
            x = "Diapause / E4.5", y = "Arrested / 8C")

    ggsave(width = 7.6, height = 7.6, filename = paste0(names(plotdat[i]), "_fc.png"),
            path = "../figs")
    ggsave(width = 7.6, height = 7.6, filename = paste0(names(plotdat[i]), "_fc.pdf"),
            path = "../figs")

}



####################################################
### for up pathways in Scognamiglio et al., 2016 ###
####################################################

plotdat <- list()
for (i in 1:length(upidx)) {
    id <- upidx[i]
    plotdat[[i]] <- cbind(dat, map_fin[[id]])
    colnames(plotdat[[i]])[3] <- "inpath"
    plotdat[[i]]$inpath <- factor(plotdat[[i]]$inpath)
    names(plotdat)[i] <- names(map_fin[id])
}


for (i in 1:length(plotdat)) {

    p <- ggplot()
    p + geom_point(data = plotdat[[i]] %>% filter(inpath == 0),
                    color = "#aab0b4", alpha = .6,
                    aes(x = mouse_fc, y = human_fc)) +
        geom_point(data = plotdat[[i]] %>% filter(inpath == 1),
                    color = "#c72a2a", alpha = .6,
                    aes(x = mouse_fc, y = human_fc)) +

        #scale_color_manual(values = c("#7C878E", "#ec4848")) +

        geom_vline(xintercept = 0, linetype = "dashed", size = 0.6) +
        geom_hline(yintercept = 0, linetype = "dashed", size = 0.6) +

        labs(title = names(plotdat[i]),
            x = "Diapause / E4.5", y = "Arrested / 8C")

    ggsave(width = 7.6, height = 7.6, filename = paste0(names(plotdat[i]), "_fc.png"),
            path = "../figs")
    ggsave(width = 7.6, height = 7.6, filename = paste0(names(plotdat[i]), "_fc.pdf"),
            path = "../figs")

}
