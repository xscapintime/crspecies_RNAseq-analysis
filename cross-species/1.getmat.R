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

human_mat <- cbind(h_arr, h_8c)

# exclude TE
human_gene_mat <- human_mat[grep("ENSG", rownames(human_mat)), ]

write.table(human_gene_mat, file = "human_readcounts.tsv",
            quote = F, sep = "\t")


# moue data
files <- list.files("../mouse_diapause/data/tecounts/", pattern = ".tsv.gz$")
m_data <- lapply(files, function(file) {
               read.table(paste0("../mouse_diapause/data/tecounts/", file),
                                sep = "\t", row.names = 1, header = F)})

spnames <- sub(".tsv.gz", "", files)
#names(m_data) <- spnames

m_mat <- Reduce(cbind, m_data)
mouse_mat <- m_mat[, c(1,3,5,7,9)]
colnames(mouse_mat) <- spnames

# exclude TE
mouse_gene_mat <- mouse_mat[grep("ENSMUSG", rownames(mouse_mat)), ]

write.table(mouse_gene_mat, file = "mouse_readcounts.tsv",
            quote = F, sep = "\t")
