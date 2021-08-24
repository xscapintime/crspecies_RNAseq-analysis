# get TSS -1000 sequence for each DEG
# ChIPseeker package
# ---------------------------------------

rm(list = ls())
library(tidyverse)


# BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene", INSTALL_opts = c('--no-lock'))
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

# options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
# BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library(BSgenome.Hsapiens.UCSC.hg38)


## ID convert
library(AnnotationDbi)
library(org.Hs.eg.db)

# type 2 up and down genes
up_t2 <- read.table("../compr-expr/common_up_deg_t2.tsv")
dn_t2 <- read.table("../compr-expr/common_dn_deg_t2.tsv")

up_t2_ez <- AnnotationDbi::select(org.Hs.eg.db, keys = up_t2$V1, columns = "ENTREZID", keytype = "SYMBOL")
dn_t2_ez <- AnnotationDbi::select(org.Hs.eg.db, keys = dn_t2$V1, columns = "ENTREZID", keytype = "SYMBOL")


## get promoter
up_trs <- transcriptsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by = "gene")[up_t2_ez$ENTREZID]
up_prom <- getPromoterSeq(up_trs, Hsapiens, upstream = 1000, downstream = 0)
names(up_prom) <- up_t2_ez$SYMBOL

dn_trs <- transcriptsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by = "gene")[dn_t2_ez$ENTREZID]
dn_prom <- getPromoterSeq(dn_trs, Hsapiens, upstream = 1000, downstream = 0)
names(dn_prom) <- dn_t2_ez$SYMBOL


## export promoter fasta
for (i in names(up_prom)) {
    writeXStringSet(up_prom[[i]], filepath = paste0("upstream_seq/common_up/", i, "_prom_1k.fasta"), format = "fasta")
}

for (i in names(dn_prom)) {
    writeXStringSet(dn_prom[[i]], filepath = paste0("upstream_seq/common_dn/", i, "_prom_1k.fasta"), format = "fasta")
}