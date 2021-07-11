# 2D MANOVA by pacakge mitch
# use logFC as input ,compare with my own result
# -------------------------

rm(list = ls())
library(tidyverse)

library(mitch)


## DESeq2 results
human_res <- read.table("../human_deg.tsv")
mouse_res <- read.table("../mouse_deg.tsv")

# need remove NA pval/qval genes
human_res <- human_res %>% na.omit()
mouse_res <- mouse_res %>% na.omit()

# NA pval in mouse data, make them be in the same length
human_res <- human_res %>% filter(row %in% mouse_res$row)



### mitch
## take 'log2FoldChange' column of each input table
any(human_res$row != mouse_res$row)

y <- cbind(human_res$log2FoldChange, mouse_res$log2FoldChange)
dimnames(y) <- list(human_res$row, c("human", "mouse"))

any(human_res$row != row.names(y))
any(mouse_res$row != row.names(y))


## all genesets
bp <- gmt_import("../gmtdata/c5.go.bp.v7.4.symbols.gmt")
mf <- gmt_import("../gmtdata/c5.go.mf.v7.4.symbols.gmt")
cc <- gmt_import("../gmtdata/c5.go.cc.v7.4.symbols.gmt")

kegg <- gmt_import("../gmtdata/c2.cp.kegg.v7.4.symbols.gmt")

# diapause gene set
load("diapuse_set.Rdata")


genesets <- c(bp, mf, cc, kegg, diapuse_set)


## mitch
res <- mitch_calc(y, genesets, priority = "significance", minsetsize = 10)
all_mitch_res <- res$enrichment_result %>% arrange(p.adjustMANOVA)



####################
#### my results ####
####################

# position score
load("position_score.Rdata")

# BH p-val
load("manova_bh_adjpval.Rdata")

my_res <- cbind(human_ps, mouse_ps, bh_adjpval) %>% data.frame %>%
                arrange(bh_adjpval)


head(my_res)
head(all_mitch_res)

# exactly the same #