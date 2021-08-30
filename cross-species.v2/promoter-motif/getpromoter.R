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

human_up_t2 <- human_res_t2 %>% dplyr::filter(log2FoldChange > 2 & padj < 0.05) %>% dplyr::select(row)
human_dn_t2 <- human_res_t2 %>% dplyr::filter(log2FoldChange < -2 & padj < 0.05) %>% dplyr::select(row)


### human arrested and mouse diapause common deg ###
### ============================================ ###
# # type 2 up and down genes
# up_t2 <- read.table("../compr-expr/common_up_deg_t2.tsv")
# dn_t2 <- read.table("../compr-expr/common_dn_deg_t2.tsv")


## match transcripts with gene symbols
## extremely slow way ##
## ==================================
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
# so the number is right, but some of them are empty
# results look weird also, dont't use them


## better way
## tirm promoter ranges
## =====================
up_prom <- promoters(src, upstream = 1000, downstream = 0, filter = SymbolFilter(human_up_t2$row), use.names=TRUE)
up_prom <- trim(up_prom)
up_prom <- split(up_prom, up_prom$symbol)

up_prom_seq <- getSeq(Hsapiens, up_prom)


dn_prom <- promoters(src, upstream = 1000, downstream = 0, filter = SymbolFilter(human_dn_t2$row), use.names=TRUE)
dn_prom <- trim(dn_prom)
dn_prom <- split(dn_prom, dn_prom$symbol)

dn_prom_seq <- getSeq(Hsapiens, dn_prom)


## export
# 624 out of 627
for (i in names(up_prom_seq)) {
    #writeXStringSet(up_prom[[i]], filepath = paste0("upstream_seq/common_up/", i, "_prom_1k.fasta"), format = "fasta")
    writeXStringSet(up_prom_seq[[i]], filepath = paste0("upstream_seq/arr_up.v2/", i, "_prom_1k.fasta"), format = "fasta")
}

# 888 out of 891
for (i in names(dn_prom_seq)) {
    #writeXStringSet(dn_prom[[i]], filepath = paste0("upstream_seq/common_dn/", i, "_prom_1k.fasta"), format = "fasta")
    writeXStringSet(dn_prom_seq[[i]], filepath = paste0("upstream_seq/arr_dn.v2/", i, "_prom_1k.fasta"), format = "fasta")
}
