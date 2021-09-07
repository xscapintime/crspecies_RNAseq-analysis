# plot motif logos
# ------------------

library(motifStack)
#library(universalmotif)

rm(list = ls())

library(tidyverse)

### common down ###
## list files
pcmfiles <- list.files("pwm/common/dn", pattern = "*.pcm")
memefiles <- list.files("pwm/common/dn", pattern = "*.meme")
homerfiles <- list.files("pwm/common/dn", pattern = "*.motif")


pcms <- lapply(pcmfiles, function(file) { t(read.table(paste0("pwm/common/dn/", file), comment.char = ">"))})
names(pcms) <- unlist(lapply(X = pcmfiles, FUN = function(x) {return(strsplit(x, split = ".pcm")[[1]][1])}))

# convert pcm to pfm
hocopfms <- lapply(pcms, pcm2pfm)


memes <- lapply(memefiles, function(file) { t(read.table(paste0("pwm/common/dn/", file), skip = 11, check.names = FALSE, comment.char = "U"))})
names(memes) <- c("CREB1_MA0018.3", "JDP2_MA0656.1", "KLF16_MA0741.1", "ATF7_MA0834.1")


homers <- lapply(homerfiles, function(file) { t(read.table(paste0("pwm/common/dn/", file), comment.char = ">"))})
names(homers) <- c("ZKSCAN1(Zf)/HepG2-ZKSCAN1-ChIP-Seq(Encode)/Homer")


motif_all <- c(hocopfms, memes, homers)

pfmlist <- list()
for (i in 1:length(motif_all)) {
    row.names(motif_all[[i]]) <- c("A","C","G","T")
    pfmlist[i] <- new("pfm", mat = as.matrix(motif_all[[i]]), name = names(motif_all[i]))
}

pdf("dn_stack.pdf") #, width = 720, height = 1500
motifStack(pfmlist, layout = "stack", ncex = 0.8, font = "Arial")
dev.off()



### common up ###
## list files
pcmfiles <- list.files("pwm/common/up", pattern = "*.pcm")
#memefiles <- list.files("pwm/common/up", pattern = "*.meme")
homerfiles <- list.files("pwm/common/up", pattern = "*.motif")


pcms <- lapply(pcmfiles, function(file) { t(read.table(paste0("pwm/common/up/", file), comment.char = ">"))})
names(pcms) <- unlist(lapply(X = pcmfiles, FUN = function(x) {return(strsplit(x, split = ".pcm")[[1]][1])}))

# convert pcm to pfm
hocopfms <- lapply(pcms, pcm2pfm)


# memes <- lapply(memefiles, function(file) { t(read.table(paste0("pwm/common/dn/", file), skip = 11, check.names = FALSE, comment.char = "U"))})
# names(memes) <- c("CREB1_MA0018.3", "JDP2_MA0656.1", "KLF16_MA0741.1", "ATF7_MA0834.1")


homers <- lapply(homerfiles, function(file) { t(read.table(paste0("pwm/common/up/", file), comment.char = ">"))})
names(homers) <- c("Nur77(NR)/K562-NR4A1-ChIP-Seq(GSE31363)/Homer")


motif_all <- c(hocopfms, homers)

pfmlist <- list()
for (i in 1:length(motif_all)) {
    row.names(motif_all[[i]]) <- c("A","C","G","T")
    pfmlist[i] <- new("pfm", mat = as.matrix(motif_all[[i]]), name = names(motif_all[i]))
}

png("t2up_stack.png", width = 500, height = 800, units = "in", res = 300)
#par(mar= c(4, 4, 2, 1))
motifStack(pfmlist, layout = "stack", font = "Arial", ncex = 1)
dev.off()
