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



###################
## plot all path ##
###################

dat <- cbind(human_res$log2FoldChange, mouse_res$log2FoldChange)
dimnames(dat) <- list(human_res$row, c("human_fc", "mouse_fc"))
dat <- data.frame(dat)


system.time({

    plotdat <- list()
    for (i in 1:length(map_fin)) {
        plotdat[[i]] <- cbind(dat, map_fin[[i]])
        colnames(plotdat[[i]])[3] <- "inpath"
        plotdat[[i]]$inpath <- factor(plotdat[[i]]$inpath)
        names(plotdat)[i] <- names(map_fin[i])
    }

})


### figure setting
library(ggplot2)
library(ggthemes)

theme_set(theme_few() + theme(legend.position = "none"))


system.time({

    for (i in 1:length(plotdat)) {

        p <- ggplot()
        p + geom_point(data = plotdat[[i]] %>% filter(inpath == 0),
                    color = "#aab0b4", alpha = .6,
                    aes(x = mouse_fc, y = human_fc)) +
            geom_point(data = plotdat[[i]] %>% filter(inpath == 1),
                    color = "#0d8168", alpha = .6,
                    aes(x = mouse_fc, y = human_fc)) +

            #scale_color_manual(values = c("#aab0b4", "#0d8168")) +

            geom_vline(xintercept = 0, linetype = "dashed", size = 0.6) +
            geom_hline(yintercept = 0, linetype = "dashed", size = 0.6) +

            labs(title = names(plotdat[i]),
                x = "Diapause / E4.5", y = "Arrested / 8C")

        ggsave(width = 7.6, height = 7.6, filename = paste0(names(plotdat[i]), "_fc.png"),
                path = "../figs/all")
    # ggsave(width = 7.6, height = 7.6, filename = paste0(names(plotdat[i]), "_fc.pdf"),
    #         path = "../figs")

    }
})
