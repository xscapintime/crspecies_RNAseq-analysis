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


# diapause sets
diaidx <- names(map_fin) %>% grep(pattern = "Diapause")



## plot data
any(human_res$row != mouse_res$row)

dat <- cbind(human_res$log2FoldChange, mouse_res$log2FoldChange)
dimnames(dat) <- list(human_res$row, c("human_fc", "mouse_fc"))
dat <- data.frame(dat)



### figure setting
library(ggplot2)
library(ggthemes)

theme_set(theme_few() + theme(legend.position = "none"))


#####################
## 4 diapause sets ##
#####################

plotdat <- list()
for (i in 1:length(diaidx)) {
    id <- diaidx[i]
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
                    color = "#EFC000FF", alpha = .6,
                    aes(x = mouse_fc, y = human_fc)) +

        #scale_color_manual(values = c("#7C878E", "#EFC000FF")) +

        geom_vline(xintercept = 0, linetype = "dashed", size = 0.6) +
        geom_hline(yintercept = 0, linetype = "dashed", size = 0.6) +

        labs(title = names(plotdat[i]),
            x = "Diapause / E4.5", y = "Arrested / 8C")

    ggsave(width = 7.6, height = 7.6, filename = paste0(names(plotdat[i]), "_fc.png"),
            path = "../figs")
    ggsave(width = 7.6, height = 7.6, filename = paste0(names(plotdat[i]), "_fc.pdf"),
            path = "../figs")

}