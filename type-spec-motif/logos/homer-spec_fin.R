# "remove" pariwise overlap motif in the common-motif list
# for homer result, only one up and dn motif are 3-types0-common
# find pairwise-overlapped motifs
# --------------------------------

rm(list = ls())

library(tidyverse)


### homer res
## list files
files_known <- list.files("../homer", pattern = "knownResults.txt", recursive = T)
files_novo <- list.files("../homer", pattern = "homerResults.txt", recursive = T)


### homer result
## homer known
homer_known <- lapply(files_known, function(file) { read.csv(paste0("../homer/", file), sep = "\t", header = T) %>% filter(P.value < 0.01)})
names(homer_known) <- unlist(lapply(X = files_known, FUN = function(x) {return(strsplit(x, split = "_MotifOutput")[[1]][1])}))


## homer de novo
homer_novo <- lapply(files_novo, function(file) { read.csv(paste0("../homer/", file), sep = "\t", header = T) %>% filter(P.value < 0.01)})
names(homer_novo) <- unlist(lapply(X = files_novo, FUN = function(x) {return(strsplit(x, split = "_MotifOutput")[[1]][1])}))


### list res
typeI_up <- union(homer_known[["typeI_up"]]$Motif.Name, homer_novo[["typeI_up"]]$Best.Match.Details)
typeII_up <- union(homer_known[["typeII_up"]]$Motif.Name, homer_novo[["typeII_up"]]$Best.Match.Details)
typeIII_up <- union(homer_known[["typeIII_up"]]$Motif.Name, homer_novo[["typeIII_up"]]$Best.Match.Details)

typeI_dn <- union(homer_known[["typeI_dn"]]$Motif.Name, homer_novo[["typeI_dn"]]$Best.Match.Details)
typeII_dn <- union(homer_known[["typeII_dn"]]$Motif.Name, homer_novo[["typeII_dn"]]$Best.Match.Details)
typeIII_dn <- union(homer_known[["typeIII_dn"]]$Motif.Name, homer_novo[["typeIII_dn"]]$Best.Match.Details)


### type I vs II
typeIvsII_up <- intersect(typeI_up, typeII_up)
typeIvsII_dn <- intersect(typeI_dn, typeII_dn)

### type I vs III
typeIvsIII_up <- intersect(typeI_up, typeIII_up)
typeIvsIII_dn <- intersect(typeI_dn, typeIII_dn)

### type II vs III
typeIIvsIII_up <- intersect(typeII_up, typeIII_up)
typeIIvsIII_dn <- intersect(typeII_dn, typeIII_dn)


### export
all_pair <- list(typeIvsII_up, typeIvsII_dn, typeIvsIII_up, typeIvsIII_dn,
                typeIIvsIII_up, typeIIvsIII_dn)
names(all_pair) <- c("typeIvsII_up", "typeIvsII_dn", "typeIvsIII_up", "typeIvsIII_dn",
                    "typeIIvsIII_up", "typeIIvsIII_dn")

for (i in 1:length(all_pair)) {
    write.table(all_pair[[i]], file = paste0("homer-", names(all_pair[i]),"_fin.txt"), quote = F, row.names = F, col.names = F, sep = "\t")
}