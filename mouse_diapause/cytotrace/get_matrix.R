## take the first column
rm(list = ls())

library(tidyverse)

setwd(paste0(getwd(), "/mouse_diapause/matrix"))

## load the tsv.gz files
files <- list.files("../data", pattern = ".tsv.gz$")
data <- lapply(files, function(file) { read.csv(paste0("../data/", file),
                                                stringsAsFactors = FALSE,
                                                sep = "\t",
                                                header = F,
                                                row.names = 1,
                                                comment.char = "#")})
data_cbind <- Reduce(cbind, data)
mat <- data_cbind[,c(1, 3, 5, 7, 9)]
colnames(mat) <- unlist(lapply(X = files, FUN = function(x) {return(strsplit(x, split = ".tsv.gz")[[1]][1])}))

write.table(mat, quote = F, sep = "\t", "diapause_e45.tsv")
save(mat, file = "mouse_diapause_te.Rdata")