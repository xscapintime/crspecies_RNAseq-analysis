rm(list = ls())

library(CytoTRACE)
library(dplyr)

setwd(paste0(getwd(), "/arrested"))

#help(package = "CytoTRACE") 

## arrested cells only
expr <- read.table("../data/norm_input.tsv", sep = "\t", row.names = 1, header = T)
genes <- row.names(expr)
stages <- colnames(expr)

expr <- as.data.frame(lapply(expr, as.numeric))
colnames(expr) <- stages
rownames(expr) <- genes

## cytotrace
#res <- CytoTRACE(expr)
#save(res, file = "cytotrace_res.Rdata")
load("../all/cytotrace_res.Rdata")

## plot res
# for E3.xx.xxxx
phe <- unlist(lapply(X = colnames(expr), FUN = function(x) {return(strsplit(x, split = ".", fixed = T)[[1]][1])}))
# for "Hs_ss_arrested_8C_rp10"
phe[grep("arrested", phe)] <- unlist(lapply(X = phe[grep("arrested", phe)], FUN = function(x) {
    return(strsplit(x, split = "_rp", fixed = F)[[1]][1])}))
# for "Zygote_E1_C1"
phe <- sub("_E[0-9]_C[0-9]*|_E[0-9]_C[0-9]*_\\w*", "", phe)

names(phe) <-  colnames(expr)
plotCytoTRACE(res, phenotype = phe)

## plot gene
plotCytoGenes(res, numOfGenes = 20)

#test <- CytoTRACE(marrow_10x_expr)
