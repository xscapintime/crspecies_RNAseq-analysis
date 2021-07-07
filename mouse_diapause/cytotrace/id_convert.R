rm(list = ls())

library(tidyverse)

setwd(paste0(getwd(), "/mouse_diapause/matrix"))

## load the matrix
load("mouse_diapause_te.Rdata")

### for id convert
library(AnnotationDbi)
library(org.Mm.eg.db)

genes <- row.names(mat)[grep(row.names(mat), pattern = "ENSMUSG")]
symb <- select(org.Mm.eg.db,
               keys = genes,
               columns = "SYMBOL",
               keytype = "ENSEMBL")

filtered <- na.omit(symb)
dup <- duplicated(filtered$SYMBOL)
filtered <- filtered[!dup, ]

## get a named matrix
idx <- row.names(mat)[grep(row.names(mat), pattern = "ENSMUSG")] %in% filtered$ENSEMBL
row.names(mat)[grep(row.names(mat), pattern = "ENSMUSG")][idx] <- filtered$SYMBOL

write.table(mat, quote = F, sep = "\t", "diapause_e45.finnamed.tsv")