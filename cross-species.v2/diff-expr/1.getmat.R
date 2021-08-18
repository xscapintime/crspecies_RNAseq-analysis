## human arrested embryo and mouse diapause embryo
## human: arrested type 2 vs morula
## human: arrested type 3 vs E4
## mouse: diapause vs E4.5

rm(list = rm())
library(tidyverse)

# human data
h_mat <- read.table("../../human_arrested/cytotrace.v2/data/no_space.tsv",
                    sep = "\t", row.names = 1, header = T)

h_morula <- h_mat %>% select(colnames(h_mat)[grep("Morula", colnames(h_mat))])
dim(h_morula)

h_e4 <- h_mat %>% select(colnames(h_mat)[grep("E4_", colnames(h_mat))])
dim(h_e4)

h_arrtype2 <- h_mat %>% select(colnames(h_mat)[grep("Arrested_[0-9]*_Type_II$", colnames(h_mat))])
dim(h_arrtype2)

h_arrtype3 <- h_mat %>% select(colnames(h_mat)[grep("Arrested_[0-9]*_Type_III$|Arrested_Z[0-9]*_Type_III$", colnames(h_mat))])
dim(h_arrtype3)


human_mat_t2 <- cbind(h_arrtype2, h_morula)
human_mat_t3 <- cbind(h_arrtype3, h_e4)


write.table(human_mat_t2, file = "human_readcounts_t2.tsv",
            quote = F, sep = "\t")
write.table(human_mat_t3, file = "human_readcounts_t3.tsv",
            quote = F, sep = "\t")