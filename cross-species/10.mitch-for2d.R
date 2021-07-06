# 2D MANOVA by pacakge mitch
# -------------------------

rm(list = ls())
library(tidyverse)

library(mitch)


## DESeq2 results
human_res <- read.table("human_deg.tsv")
mouse_res <- read.table("mouse_deg.tsv")

# need remove NA pval/qval genes
human_res <- human_res %>% na.omit()
mouse_res <- mouse_res %>% na.omit()

# NA pval in mouse data, make them be in the same length
human_res <- human_res %>% filter(row %in% mouse_res$row)



### mitch
x <- list("human" = human_res, "mouse" = mouse_res)

y <- mitch_import(x, DEtype = "DESeq2", geneIDcol = "row")


## GO BP
genesets <- gmt_import("gmtdata/c5.go.bp.v7.4.symbols.gmt")


res <- mitch_calc(y, genesets, priority = "significance", minsetsize = 10)
mitch_bp <- res$enrichment_result

write.table(mitch_bp, file = "mitch_res_bp.txt", quote = F, sep = "\t")

#mitch_plots(res, outfile = "mitch_charts.pdf")


## GO MF
genesets <- gmt_import("gmtdata/c5.go.mf.v7.4.symbols.gmt")

res <- mitch_calc(y, genesets, priority = "significance", minsetsize = 10)
mitch_mf <- res$enrichment_result

write.table(mitch_mf, file = "mitch_res_mf.txt", quote = F, sep = "\t")


## GO CC
genesets <- gmt_import("gmtdata/c5.go.cc.v7.4.symbols.gmt")

res <- mitch_calc(y, genesets, priority = "significance", minsetsize = 10)
mitch_cc <- res$enrichment_result

write.table(mitch_cc, file = "mitch_res_cc.txt", quote = F, sep = "\t")
