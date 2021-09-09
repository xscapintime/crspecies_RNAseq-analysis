# for meme results
# find pairwise-overlapped motifs
# ---------------------------------------------------------------

rm(list = ls())

library(tidyverse)

## load tables
typeI_up <- read.csv("meme-typeI_spec_up.txt", sep = " ", header = F)
typeII_up <- read.csv("meme-typeII_spec_up.txt", sep = " ", header = F)
typeIII_up <- read.csv("meme-typeIII_spec_up.txt", sep = " ", header = F)

typeI_dn <- read.csv("meme-typeI_spec_dn.txt", sep = " ", header = F)
typeII_dn <- read.csv("meme-typeII_spec_dn.txt", sep = " ", header = F)
typeIII_dn <- read.csv("meme-typeIII_spec_dn.txt", sep = " ", header = F)


### pairwise-overlapped motifs
### type I vs II
typeIvsII_up <- inner_join(typeI_up, typeII_up, by = c("V1" = "V1"), suffix = c(".I", ".II")) %>% arrange(V2.I, V2.II)
typeIvsII_dn <- inner_join(typeI_dn, typeII_dn, by = c("V1" = "V1"), suffix = c(".I", ".II")) %>% arrange(V2.I, V2.II)

### type I vs III
typeIvsIII_up <- inner_join(typeI_up, typeIII_up, by = c("V1" = "V1"), suffix = c(".I", ".III")) %>% arrange(V2.I, V2.III)
typeIvsIII_dn <- inner_join(typeI_dn, typeIII_dn, by = c("V1" = "V1"), suffix = c(".I", ".III")) %>% arrange(V2.I, V2.III)

### type II vs III
typeIIvsIII_up <- inner_join(typeII_up, typeIII_up, by = c("V1" = "V1"), suffix = c(".II", ".III")) %>% arrange(V2.II, V2.III)
typeIIvsIII_dn <- inner_join(typeII_up, typeIII_up, by = c("V1" = "V1"), suffix = c(".II", ".III")) %>% arrange(V2.II, V2.III)


### export
all_pair <- list(typeIvsII_up, typeIvsII_dn, typeIvsIII_up, typeIvsIII_dn,
                typeIIvsIII_up, typeIIvsIII_dn)
names(all_pair) <- c("typeIvsII_up", "typeIvsII_dn", "typeIvsIII_up", "typeIvsIII_dn",
                    "typeIIvsIII_up", "typeIIvsIII_dn")

for (i in 1:length(all_pair)) {
    write.table(all_pair[[i]], file = paste0("meme-", names(all_pair[i]),"_fin.txt"), quote = F, row.names = F, col.names = F, sep = "\t")
}