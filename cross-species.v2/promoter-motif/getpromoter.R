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


## ID convert
library(AnnotationDbi)
library(org.Hs.eg.db)

up_t2_ez <- AnnotationDbi::select(src, keys = human_up_t2$row, columns = "entrez", keytype = "symbol") %>% distinct()
dn_t2_ez <- AnnotationDbi::select(TxDb.Hsapiens.UCSC.hg38.knownGene, keys = human_dn_t2$row, columns = "TXNAME", keytype = "TXNAME")


library(biomaRt)
human <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
human_go_id <- getBM(attributes = c("hgnc_symbol", "entrezgene_id"),
                        filters = "hgnc_symbol",
                        values = human_up_t2$row, mart = human)



## get promoter
all_trs <- transcriptsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by = "gene")#[up_t2_ez$entrez]

up_prom <- getPromoterSeq(up_trs, Hsapiens, upstream = 1000, downstream = 0)
names(up_prom) <- up_t2_ez$SYMBOL

dn_trs <- transcriptsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by = "gene")[dn_t2_ez$ENTREZID]
dn_prom <- getPromoterSeq(dn_trs, Hsapiens, upstream = 1000, downstream = 0)
names(dn_prom) <- dn_t2_ez$SYMBOL





up_trs <- GRangesList()
for(i in 1:nrow(human_up_t2)) { 
    up_trs[[i]] <- transcripts(src, filter = SymbolFilter(human_dn_t2$row[i]))
    up_trs[[i]] <- up_trs[[i]][!seqnames(up_trs[[i]]) %>% stringr::str_detect("_")]
}
names(up_trs) <- up_t2_ez$symbol

up_prom <- getPromoterSeq(up_trs, Hsapiens, upstream = 1000)



## export promoter fasta
for (i in names(up_prom)) {
    writeXStringSet(up_prom[[i]], filepath = paste0("upstream_seq/common_up/", i, "_prom_1k.fasta"), format = "fasta")
}

for (i in names(dn_prom)) {
    writeXStringSet(dn_prom[[i]], filepath = paste0("upstream_seq/common_dn/", i, "_prom_1k.fasta"), format = "fasta")
}