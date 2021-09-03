BiocManager::install("motifStack", INSTALL_opts = c('--no-lock'))

library(motifStack)

pcm <- read.table(file.path(find.package("motifStack"),
                            "extdata", "bin_SOLEXA.pcm"))
pcm <- pcm[,3:ncol(pcm)]
rownames(pcm) <- c("A","C","G","T")
motif <- new("pcm", mat=as.matrix(pcm), name="bin_SOLEXA")


motif <- new("pcm", mat=as.matrix(pcm), name="bin_SOLEXA")


browseVignettes("motifStack")
