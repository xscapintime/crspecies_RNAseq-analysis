# remove pariwise overlap motif
# for meme results
# --------------------------------

rm(list = ls())

library(tidyverse)


## load tables
typeI_up <- read.csv("meme-typeI_spec_up.txt", sep = " ", header = F)
typeII_up <- read.csv("meme-typeII_spec_up.txt", sep = " ", header = F)
typeIII_up <- read.csv("meme-typeIII_spec_up.txt", sep = " ", header = F)

typeI_dn <- read.csv("meme-typeI_spec_dn.txt", sep = " ", header = F)
typeII_dn <- read.csv("meme-typeII_spec_dn.txt", sep = " ", header = F)
typeIII_dn <- read.csv("meme-typeIII_spec_dn.txt", sep = " ", header = F)


## pariwise excluding
todel_up <- c(intersect(typeI_up$V1, typeII_up$V1), intersect(typeI_up$V1, typeIII_up$V1), intersect(typeIII_up$V1, typeII_up$V1))
todel_dn <- c(intersect(typeI_dn$V1, typeII_dn$V1), intersect(typeI_dn$V1, typeIII_dn$V1), intersect(typeIII_dn$V1, typeII_dn$V1))


## remove
typeI_up <- typeI_up %>% filter(!V1 %in% todel_up) %>% arrange(V2)
typeII_up <- typeII_up %>% filter(!V1 %in% todel_up) %>% arrange(V2)
typeIII_up <- typeIII_up %>% filter(!V1 %in% todel_up) %>% arrange(V2)

typeI_dn <- typeI_dn %>% filter(!V1 %in% todel_dn) %>% arrange(V2)
typeII_dn <- typeII_dn %>% filter(!V1 %in% todel_dn) %>% arrange(V2)
typeIII_dn <- typeIII_dn %>% filter(!V1 %in% todel_dn) %>% arrange(V2)


allspec <- list(typeI_up, typeII_up, typeIII_up, typeI_dn, typeII_dn, typeIII_dn)
names(allspec) <- c("typeI_up", "typeII_up", "typeIII_up", "typeI_dn", "typeII_dn", "typeIII_dn")

for (i in 1:length(allspec)) {
    write.table(allspec[[i]], file = paste0("meme-", names(allspec[i]),"_fin.txt"), quote = F, row.names = F, col.names = F, sep = "\t")
}

