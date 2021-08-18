rm(list = ls())

library(CytoTRACE)
library(tidyverse)

#setwd(paste0(getwd(), "/arrested"))

#help(package = "CytoTRACE") 

## arrested cells only
expr <- read.table("data/no_space.tsv", sep = "\t", row.names = 1, header = T)
genes <- row.names(expr)
stages <- colnames(expr)

expr <- as.data.frame(lapply(expr, as.numeric))
colnames(expr) <- stages
rownames(expr) <- genes

## cytotrace
res <- CytoTRACE(expr)
save(res, file = "cytotrace_res.Rdata")

#########======================#####################
load("cytotrace_res.Rdata")

## plot res
# group the phenotypes

phe <- colnames(expr)
phe[grep("Type", phe)] <- unlist(lapply(X = phe[grep("Type", phe)], FUN = function(x) {
    return(strsplit(x, split = "_Type_", fixed = F)[[1]][2])}))

#phe <- colnames(expr)[grep("Arrested", colnames(expr), invert = T)]

phe <- ifelse(grepl("_", phe) == TRUE,
            unlist(lapply(X = phe, FUN = function(x) {return(strsplit(x, split = "_", fixed = T)[[1]][1])})),
            phe)

phe[grep("Arrested", colnames(expr))] <- paste0("Arrested_Type_", phe[grep("Arrested", colnames(expr))])
phe[grep("Late", phe)] <- paste0(phe[grep("Late", phe)], "_blastocyst")

# 2C, 4C, 8C
phe <- sub("X", "", phe)

names(phe) <- colnames(expr)
plotCytoTRACE(res, phenotype = phe)


## plot gene
plotCytoGenes(res, numOfGenes = 20)

#test <- CytoTRACE(marrow_10x_expr)
