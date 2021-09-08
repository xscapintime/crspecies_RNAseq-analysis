# plot motif logos
# all type-specific ones
# ------------------

library(motifStack)
#library(universalmotif)

rm(list = ls())

library(tidyverse)

###########################################
###=============== typeI ===============###
###########################################
## up
pcmfiles <- list.files("pwm/typeI/up", pattern = "*.pcm")
pfmfiles <- list.files("pwm/typeI/up", pattern = "*.pfm")
memefiles <- list.files("pwm/typeI/up", pattern = "*.meme")

# read
# hoco
pcms <- lapply(pcmfiles, function(file) { t(read.table(paste0("pwm/typeI/up/", file), comment.char = ">"))})
names(pcms) <- unlist(lapply(X = pcmfiles, FUN = function(x) {return(strsplit(x, split = ".pcm")[[1]][1])}))

# jaspar pfm looks like pcm...
pfms <- lapply(pfmfiles, function(file) { read.table(paste0("pwm/typeI/up/", file), comment.char = ">")})
names(pfms) <- unlist(lapply(X = pfmfiles, FUN = function(x) {return(strsplit(x, split = ".pfm")[[1]][1])}))

# from jolma2013
memes <- lapply(memefiles, function(file) { t(read.table(paste0("pwm/typeI/up/", file), skip = 3, check.names = FALSE))})
names(memes) <- unlist(lapply(X = memefiles, FUN = function(x) {return(strsplit(x, split = ".meme")[[1]][1])}))


# convert pcm to pfm
hocopfms <- lapply(pcms, pcm2pfm)
jsppfms <- lapply(pfms, pcm2pfm)


# plot
motif_all <- c(hocopfms, jsppfms, memes)

pfmlist <- list()
for (i in 1:length(motif_all)) {
    row.names(motif_all[[i]]) <- c("A","C","G","T")
    pfmlist[i] <- new("pfm", mat = as.matrix(motif_all[[i]]), name = names(motif_all[i]))
}
names(pfmlist) <- names(motif_all)

#pdf("typeI_up_stack_full.pdf", width = 8, height = 28)
pdf("typeI_up_tree_full.pdf", width = 10, height = 28)
#par(mar = c(4.1, 4.1, 2.1, 4.1))
motifStack(pfmlist, layout = "tree", ncex = 0.8, font = "Arial", reorder = T)
while (!is.null(dev.list()))  dev.off()


pfmlist_drop <- pfmlist[c("ELF1_ETS_2", "ZNF24_MA1124.1", "SP4_MA0685.1", "KLF9_MA1107.1",
                            "KLF14_MA0740.1", "ELF4_MA0641.1", "ELF3_MA0640.1", "ZNF76_HUMAN.H11MO.0.C",
                            "ZN121_HUMAN.H11MO.0.C", "TYY1_HUMAN.H11MO.0.A", "THA11_HUMAN.H11MO.0.B",
                            "MYCN_HUMAN.H11MO.0.A", "KLF4_HUMAN.H11MO.0.A", "ELK1_HUMAN.H11MO.0.B")]

pdf("typeI_up_tree_drop.pdf", width = 10, height = 18)
motifStack(pfmlist_drop, layout = "tree", ncex = 0.8, font = "Arial", reorder = T)
while (!is.null(dev.list()))  dev.off()


## dn
pcmfiles <- list.files("pwm/typeI/dn", pattern = "*.pcm")
pfmfiles <- list.files("pwm/typeI/dn", pattern = "*.pfm")
memefiles <- list.files("pwm/typeI/dn", pattern = "*.meme")

# read
# hoco
pcms <- lapply(pcmfiles, function(file) { t(read.table(paste0("pwm/typeI/dn/", file), comment.char = ">"))})
names(pcms) <- unlist(lapply(X = pcmfiles, FUN = function(x) {return(strsplit(x, split = ".pcm")[[1]][1])}))

# jaspar pfm looks like pcm...
pfms <- lapply(pfmfiles, function(file) { read.table(paste0("pwm/typeI/dn/", file), comment.char = ">")})
names(pfms) <- unlist(lapply(X = pfmfiles, FUN = function(x) {return(strsplit(x, split = ".pfm")[[1]][1])}))

# from jolma2013
memes <- lapply(memefiles, function(file) { t(read.table(paste0("pwm/typeI/dn/", file), skip = 3, check.names = FALSE))})
names(memes) <- unlist(lapply(X = memefiles, FUN = function(x) {return(strsplit(x, split = ".meme")[[1]][1])}))


# convert pcm to pfm
hocopfms <- lapply(pcms, pcm2pfm)
jsppfms <- lapply(pfms, pcm2pfm)


# plot
motif_all <- c(hocopfms, jsppfms, memes)

pfmlist <- list()
for (i in 1:length(motif_all)) {
    row.names(motif_all[[i]]) <- c("A","C","G","T")
    pfmlist[i] <- new("pfm", mat = as.matrix(motif_all[[i]]), name = names(motif_all[i]))
}
names(pfmlist) <- names(motif_all)

#pdf("typeI_up_stack_full.pdf", width = 8, height = 28)
pdf("typeI_dn_tree_full.pdf", width = 10, height = 11)
#par(mar = c(4.1, 4.1, 2.1, 4.1))
motifStack(pfmlist, layout = "tree", ncex = 0.8, font = "Arial", reorder = T)
while (!is.null(dev.list()))  dev.off()



###########################################
###=============== typeII ==============###
###########################################
## up
pcmfiles <- list.files("pwm/typeII/up", pattern = "*.pcm")
pfmfiles <- list.files("pwm/typeII/up", pattern = "*.pfm")
# memefiles <- list.files("pwm/typeII/up", pattern = "*.meme")

# read
# hoco
pcms <- lapply(pcmfiles, function(file) { t(read.table(paste0("pwm/typeII/up/", file), comment.char = ">"))})
names(pcms) <- unlist(lapply(X = pcmfiles, FUN = function(x) {return(strsplit(x, split = ".pcm")[[1]][1])}))

# jaspar pfm looks like pcm...
pfms <- lapply(pfmfiles, function(file) { read.table(paste0("pwm/typeII/up/", file), comment.char = ">")})
names(pfms) <- unlist(lapply(X = pfmfiles, FUN = function(x) {return(strsplit(x, split = ".pfm")[[1]][1])}))

# from jolma2013
# memes <- lapply(memefiles, function(file) { t(read.table(paste0("pwm/typeII/up/", file), skip = 3, check.names = FALSE))})
# names(memes) <- unlist(lapply(X = memefiles, FUN = function(x) {return(strsplit(x, split = ".meme")[[1]][1])}))


# convert pcm to pfm
hocopfms <- lapply(pcms, pcm2pfm)
jsppfms <- lapply(pfms, pcm2pfm)


# plot
motif_all <- c(hocopfms, jsppfms)

pfmlist <- list()
for (i in 1:length(motif_all)) {
    row.names(motif_all[[i]]) <- c("A","C","G","T")
    pfmlist[i] <- new("pfm", mat = as.matrix(motif_all[[i]]), name = names(motif_all[i]))
}
names(pfmlist) <- names(motif_all)

#pdf("typeI_up_stack_full.pdf", width = 8, height = 28)
pdf("typeII_up_tree_full.pdf", width = 10, height = 9)
#par(mar = c(4.1, 4.1, 2.1, 4.1))
motifStack(pfmlist, layout = "tree", ncex = 0.8, font = "Arial", reorder = T)
while (!is.null(dev.list()))  dev.off()



## dn
pcmfiles <- list.files("pwm/typeII/dn", pattern = "*.pcm")
pfmfiles <- list.files("pwm/typeII/dn", pattern = "*.pfm")
memefiles <- list.files("pwm/typeII/dn", pattern = "*.meme")

# read
# hoco
pcms <- lapply(pcmfiles, function(file) { t(read.table(paste0("pwm/typeII/dn/", file), comment.char = ">"))})
names(pcms) <- unlist(lapply(X = pcmfiles, FUN = function(x) {return(strsplit(x, split = ".pcm")[[1]][1])}))

# jaspar pfm looks like pcm...
pfms <- lapply(pfmfiles, function(file) { read.table(paste0("pwm/typeII/dn/", file), comment.char = ">")})
names(pfms) <- unlist(lapply(X = pfmfiles, FUN = function(x) {return(strsplit(x, split = ".pfm")[[1]][1])}))

# from jolma2013
memes <- lapply(memefiles, function(file) { t(read.table(paste0("pwm/typeII/dn/", file), skip = 3, check.names = FALSE))})
names(memes) <- unlist(lapply(X = memefiles, FUN = function(x) {return(strsplit(x, split = ".meme")[[1]][1])}))


# convert pcm to pfm
hocopfms <- lapply(pcms, pcm2pfm)
jsppfms <- lapply(pfms, pcm2pfm)


# plot
motif_all <- c(hocopfms, jsppfms, memes)

pfmlist <- list()
for (i in 1:length(motif_all)) {
    row.names(motif_all[[i]]) <- c("A","C","G","T")
    pfmlist[i] <- new("pfm", mat = as.matrix(motif_all[[i]]), name = names(motif_all[i]))
}
names(pfmlist) <- names(motif_all)

#pdf("typeI_up_stack_full.pdf", width = 8, height = 28)
pdf("typeII_dn_tree_full.pdf", width = 10, height = 27)
#par(mar = c(4.1, 4.1, 2.1, 4.1))
motifStack(pfmlist, layout = "tree", ncex = 0.8, font = "Arial", reorder = T)
while (!is.null(dev.list()))  dev.off()




###########################################
###============== typeIII ==============###
###########################################
## up
pcmfiles <- list.files("pwm/typeIII/up", pattern = "*.pcm")
# pfmfiles <- list.files("pwm/typeIII/up", pattern = "*.pfm")
# memefiles <- list.files("pwm/typeIII/up", pattern = "*.meme")

# read
# hoco
pcms <- lapply(pcmfiles, function(file) { t(read.table(paste0("pwm/typeIII/up/", file), comment.char = ">"))})
names(pcms) <- unlist(lapply(X = pcmfiles, FUN = function(x) {return(strsplit(x, split = ".pcm")[[1]][1])}))

# jaspar pfm looks like pcm...
# pfms <- lapply(pfmfiles, function(file) { read.table(paste0("pwm/typeIII/up/", file), comment.char = ">")})
# names(pfms) <- unlist(lapply(X = pfmfiles, FUN = function(x) {return(strsplit(x, split = ".pfm")[[1]][1])}))

# from jolma2013
# memes <- lapply(memefiles, function(file) { t(read.table(paste0("pwm/typeIII/up/", file), skip = 3, check.names = FALSE))})
# names(memes) <- unlist(lapply(X = memefiles, FUN = function(x) {return(strsplit(x, split = ".meme")[[1]][1])}))


# convert pcm to pfm
hocopfms <- lapply(pcms, pcm2pfm)
# jsppfms <- lapply(pfms, pcm2pfm)


# plot
motif_all <- c(hocopfms)

pfmlist <- list()
for (i in 1:length(motif_all)) {
    row.names(motif_all[[i]]) <- c("A","C","G","T")
    pfmlist[i] <- new("pfm", mat = as.matrix(motif_all[[i]]), name = names(motif_all[i]))
}
names(pfmlist) <- names(motif_all)

#pdf("typeIII_up_stack_full.pdf", width = 8, height = 28)
pdf("typeIII_up_full.pdf", width = 8, height = 4) # there's only one
# motifStack(pfmlist, layout = "stack", ncex = 0.8, font = "Arial", reorder = T)
plot(pfmlist[["ZSC22_HUMAN.H11MO.0.C"]], font = "Arial")
while (!is.null(dev.list()))  dev.off()



## dn
# pcmfiles <- list.files("pwm/typeIII/dn", pattern = "*.pcm")
pfmfiles <- list.files("pwm/typeIII/dn", pattern = "*.pfm")
# memefiles <- list.files("pwm/typeIII/dn", pattern = "*.meme")

# read
# hoco
# pcms <- lapply(pcmfiles, function(file) { t(read.table(paste0("pwm/typeIII/dn/", file), comment.char = ">"))})
# names(pcms) <- unlist(lapply(X = pcmfiles, FUN = function(x) {return(strsplit(x, split = ".pcm")[[1]][1])}))

# jaspar pfm looks like pcm...
pfms <- lapply(pfmfiles, function(file) { read.table(paste0("pwm/typeIII/dn/", file), comment.char = ">")})
names(pfms) <- unlist(lapply(X = pfmfiles, FUN = function(x) {return(strsplit(x, split = ".pfm")[[1]][1])}))

# from jolma2013
# memes <- lapply(memefiles, function(file) { t(read.table(paste0("pwm/typeIII/dn/", file), skip = 3, check.names = FALSE))})
# names(memes) <- unlist(lapply(X = memefiles, FUN = function(x) {return(strsplit(x, split = ".meme")[[1]][1])}))


# convert pcm to pfm
# hocopfms <- lapply(pcms, pcm2pfm)
jsppfms <- lapply(pfms, pcm2pfm)


# plot
motif_all <- c(jsppfms)

pfmlist <- list()
for (i in 1:length(motif_all)) {
    row.names(motif_all[[i]]) <- c("A","C","G","T")
    pfmlist[i] <- new("pfm", mat = as.matrix(motif_all[[i]]), name = names(motif_all[i]))
}
names(pfmlist) <- names(motif_all)

# png("test.png")
pdf("typeIII_dn_stack_full.pdf", width = 8, height = 6)
# pdf("test.pdf")
# pdf("typeIII_dn_tree_full.pdf", width = 10, height = 8)
plotMotifLogoStack(pfmlist, font = "Arial", ncex = 0.8) # font arg of motifStack always down

while (!is.null(dev.list()))  dev.off()

