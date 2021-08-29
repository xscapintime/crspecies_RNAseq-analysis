# get TSS -1000 sequence for each DEG
# try ensembledb package
# -----------------------------------------

options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
BiocManager::install("ensembldb")
BiocManager::install("EnsDb.Hsapiens.v86")

library(ensembldb)
library(EnsDb.Hsapiens.v86)

edb <- EnsDb.Hsapiens.v86

transcripts(edb, filter = GeneNameFilter("BCL2L11"))

transcriptsBy(edb, filter = GeneNameFilter("BCL2L11"))