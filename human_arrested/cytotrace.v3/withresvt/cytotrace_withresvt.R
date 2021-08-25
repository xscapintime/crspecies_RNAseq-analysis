rm(list = ls())

library(CytoTRACE)
library(tidyverse)

#setwd(paste0(getwd(), "/arrested"))

#help(package = "CytoTRACE")

## arrested cells only
expr <- read.table("../data/norm_input-withresveratrol.tsv/norm_input.tsv", sep = "\t", row.names = 1, header = T)
genes <- row.names(expr)
stages <- colnames(expr)

expr <- as.data.frame(lapply(expr, as.numeric))
colnames(expr) <- sub("X", "", stages)
rownames(expr) <- genes

## cytotrace
res <- CytoTRACE(expr)
save(res, file = "cytotrace_res.Rdata")

#########======================#####################
load("cytotrace_res.Rdata")

## meta data for group
meta <- read.csv("../data/norm_input-withresveratrol.tsv/sample_metadata.tsv", header = T, sep = "\t")

# check names
setdiff(meta$original_name, colnames(res$exprMatrix))
idx <- match(meta$original_name, colnames(res$exprMatrix))

phe <- meta$phenotype
names(phe) <- colnames(res$exprMatrix)[idx]

phe <- sub("\\(Untreated\\) ", "", phe)


## plot res
## cytotrace plot
plotCytoTRACE(res, phenotype = phe)


## plot gene
plotCytoGenes(res, numOfGenes = 20)


# group the phenotypes
# string stuff
# useless when got the meta table
# -----------------------
phe <- colnames(expr)

# 1) for E3, E4, E5, E6, E7, split by '.'
phe <- unlist(lapply(X = colnames(expr), FUN = function(x) {return(strsplit(x, split = ".", fixed = T)[[1]][1])}))

# 2a) for arrested with resvt
phe[grep("_res_", phe)] <- unlist(lapply(X = phe[grep("_res_", phe)], FUN = function(x) {
    return(strsplit(x, split = "_rp", fixed = F)[[1]][1])}))

# 2b) for arrested 8C, what's ZZU?
phe[grep("_8C_", phe)] <- unlist(lapply(X = phe[grep("_8C_", phe)], FUN = function(x) {
    return(strsplit(x, split = "_rp", fixed = F)[[1]][1])}))
# for _ss_
phe[grep("_8C$", phe)] <- "Hs_8C_arrested"

# 3) for Late_blastocyst
phe[grep("Late", phe)] <-  unlist(lapply(X = phe[grep("Late", phe)], FUN = function(x) {
    return(strsplit(x, split = "_E", fixed = F)[[1]][1])}))

# 4) for Epiblast
phe[grep("Epiblast", phe)] <- "Epiblast"

# 5) for Trophectoderm
phe[grep("Trophectoderm", phe)] <- "Trophectoderm"

# 6) for Endoderm
phe[grep("Endoderm", phe)] <- "Endoderm"

# 7) for Morula
phe[grep("Morula", phe)] <- "Morula"

# 8) for Oocyte
phe[grep("Oocyte", phe)] <- "Oocyte"

# 9) for Morula
phe[grep("Zygote", phe)] <- "Zygote"

# 10) for 2C, 4C, 8C
phe[grep("X", phe)] <- sub("X", "", unlist(lapply(X = phe[grep("X", phe)], FUN = function(x) {
    return(strsplit(x, split = "_E", fixed = F)[[1]][1])})))


# complete grouping!
names(phe) <- colnames(expr)




#test <- CytoTRACE(marrow_10x_expr)
