# get TSS -1000 sequence for each DEG
# human arrested embryo type 2 DEG
# ---------------------------------------

rm(list = ls())
library(tidyverse)


# BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene", INSTALL_opts = c('--no-lock'))
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# BiocManager::install("Organism.dplyr")
library(Organism.dplyr)
# src <- src_organism("TxDb.Hsapiens.UCSC.hg38.knownGene")
## create an src_organism object using the most recent TxDb version available
src <- src_ucsc("Homo sapiens")


# options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
# BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library(BSgenome.Hsapiens.UCSC.hg38)


## human arrested type DEG table
human_res_t2 <- read.table("../diff-expr/human_deg_t2.tsv")

human_up_t2 <- human_res_t2 %>% filter(log2FoldChange > 2 & padj < 0.05) %>% dplyr::select(row)
human_dn_t2 <- human_res_t2 %>% filter(log2FoldChange < -2 & padj < 0.05) %>% dplyr::select(row)


### human arrested and mouse diapause common deg ###
### ============================================ ###
# # type 2 up and down genes
# up_t2 <- read.table("../compr-expr/common_up_deg_t2.tsv")
# dn_t2 <- read.table("../compr-expr/common_dn_deg_t2.tsv")


## match transcripts with gene symbols

system.time({
    up_trs <- GRangesList()
    for(i in 1:nrow(human_up_t2)) {
        up_trs[[i]] <- transcripts(src, filter = SymbolFilter(human_up_t2$row[i]))
        up_trs[[i]] <- up_trs[[i]][!seqnames(up_trs[[i]]) %>% stringr::str_detect("_")]
    }
})
names(up_trs) <- human_up_t2$row

up_prom <- getPromoterSeq(up_trs, Hsapiens, upstream = 1000)


system.time({
    dn_trs <- GRangesList()
    for(i in 1:nrow(human_dn_t2)) {
        dn_trs[[i]] <- transcripts(src, filter = SymbolFilter(human_dn_t2$row[i]))
        dn_trs[[i]] <- dn_trs[[i]][!seqnames(dn_trs[[i]]) %>% stringr::str_detect("_")]
    }
})
names(dn_trs) <- human_dn_t2$row

dn_prom <- getPromoterSeq(dn_trs, Hsapiens, upstream = 1000)


## export promoter fasta
for (i in names(up_prom)) {
    #writeXStringSet(up_prom[[i]], filepath = paste0("upstream_seq/common_up/", i, "_prom_1k.fasta"), format = "fasta")
    writeXStringSet(up_prom[[i]], filepath = paste0("upstream_seq/arr_up/", i, "_prom_1k.fasta"), format = "fasta")
}

for (i in names(dn_prom)) {
    #writeXStringSet(dn_prom[[i]], filepath = paste0("upstream_seq/common_dn/", i, "_prom_1k.fasta"), format = "fasta")
    writeXStringSet(dn_prom[[i]], filepath = paste0("upstream_seq/arr_dn/", i, "_prom_1k.fasta"), format = "fasta")
}



tmp <- GRangesList()
    for(i in 1:nrow(human_dn_t2)) {
        dn_trs[[i]] <- transcripts(src, filter = SymbolFilter(human_dn_t2$row[i]))
        #dn_trs[[i]] <- dn_trs[[i]][!seqnames(dn_trs[[i]]) %>% stringr::str_detect("_")]
    }


