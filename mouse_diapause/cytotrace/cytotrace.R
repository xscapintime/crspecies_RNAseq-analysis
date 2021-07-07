rm(list = ls())

library(CytoTRACE)
library(tidyverse)

setwd(paste0(getwd(), "/mouse_diapause/matrix"))

help(package = "CytoTRACE")

expr <- read.table("diapause_e45.finnamed.tsv", sep = "\t", row.names = 1, header = T)
genes <- row.names(expr)
stages <- colnames(expr)

expr <- as.data.frame(lapply(expr, as.numeric))
colnames(expr) <- stages
rownames(expr) <- genes

## cytotrace
res <- CytoTRACE(expr)
#save(res, file = "cytotrace_res.Rdata")

## plot res
## merge reps
phe <- unlist(lapply(X = colnames(expr), FUN = function(x) {return(strsplit(x, split = "_R", fixed = T)[[1]][1])}))
names(phe) <-  colnames(expr)

#error: perplexity is too large for the number of samples
### when sample < 100, perplexity = 1
plotCytoTRACE_mod <- plotCytoTRACE
trace("plotCytoTRACE_mod", edit = TRUE)

plotCytoTRACE_mod(res, phenotype = phe)

## seprate
phe_r <- colnames(expr)
names(phe_r) <-  colnames(expr)

plotCytoTRACE_mod(res, phenotype = phe_r, outputDir = "./reps/")


## cytogenes
plotCytoGenes(res, numOfGenes = 20)
