## human arrested embryo and mouse diapause embryo
## human: arrested vs 8C
## mouse: diapause vs E4.5

rm(list = rm())
library(tidyverse)

# human data
h_mat <- read.table("../human_arrested/data/norm_input.tsv",
                    sep = "\t", row.names = 1, header = T)
h_8c <- h_mat %>% select(colnames(h_mat)[grep("X8C_E[0-9]_C[0-9]", colnames(h_mat))])
dim(h_8c)

h_arr <- h_mat %>% select(colnames(h_mat)[grep("Hs_ss_arrested_8C_*", colnames(h_mat))])
dim(h_arr)

# moue data
files <- list.files("../mouse_diapause/data/tecounts/", pattern = ".tsv.gz$")
m_data <- lapply(files, function(file) {
               read.table(paste0("../mouse_diapause/data/tecounts/", file),
                                sep = "\t", row.names = 1, header = F)[, 1]})

m_mat <- Reduce(cbind, m_data)