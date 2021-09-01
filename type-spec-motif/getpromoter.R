# get TSS -1000 sequence for each DEG
# for all 3 types
# ---------------------------------------

rm(list = ls())
library(tidyverse)

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(Organism.dplyr)
library(BSgenome.Hsapiens.UCSC.hg38)

src <- src_ucsc("Homo sapiens")


## human arrested type DEG table
# type 1
up_t1 <- read.csv("deg/TypeI_DE/deseq.typeIvs2C.up.tsv", header = T, sep = "\t", row.names = 1)
dn_t1 <- read.csv("deg/TypeI_DE/deseq.typeIvs2C.dn.tsv", header = T, sep = "\t", row.names = 1)

# type 2
up_t2 <- read.csv("deg/TypeII_DE/deseq.TypeIIvsMorula.up.tsv", header = T, sep = "\t", row.names = 1)
dn_t2 <- read.csv("deg/TypeII_DE/deseq.TypeIIvsMorula.dn.tsv", header = T, sep = "\t", row.names = 1)

# type 3
up_t3 <- read.csv("deg/TypeIII_DE/deseq.TypeIIIvsE4.up.tsv", header = T, sep = "\t", row.names = 1)
dn_t3 <- read.csv("deg/TypeIII_DE/deseq.TypeIIIvsE4.dn.tsv", header = T, sep = "\t", row.names = 1)



## gene symbol list
up_t1_syb <- up_t1[grep("ENSG", row.names(up_t1)),]$name
dn_t1_syb <- dn_t1[grep("ENSG", row.names(dn_t1)),]$name

up_t2_syb <- up_t2[grep("ENSG", row.names(up_t2)),]$name
dn_t2_syb <- dn_t2[grep("ENSG", row.names(dn_t2)),]$name

up_t3_syb <- up_t3[grep("ENSG", row.names(up_t3)),]$name
dn_t3_syb <- dn_t3[grep("ENSG", row.names(dn_t3)),]$name


## quick promoter retreiving
source("rush_promo.R")


genelist <- list(up_t1_syb, up_t2_syb, up_t3_syb, dn_t1_syb, dn_t2_syb, dn_t3_syb)
names(genelist) <- c("typeI_up", "typeII_up", "typeIII_up", "typeI_dn", "typeII_dn", "typeIII_dn")


## get promoter fasta
for (i in 1:6) {
    prom_seq <- rush_promo(db = src, upstream = 1000, downstream = 0, genes = genelist[[i]], sequence = T)
    writeXStringSet(unlist(prom_seq), filepath = paste0("prom_fatsa/", names(genelist)[i], "_prom_1k.fasta"), format = "fasta")
}


## get promoter bed
for (i in 1:6) {
    prom <- rush_promo(db = src, upstream = 1000, downstream = 0, genes = genelist[[i]], sequence = F)
    rtracklayer::export.bed(unlist(prom), paste0("prom_bed/", names(genelist)[i], "_prom_1k.bed"))
}