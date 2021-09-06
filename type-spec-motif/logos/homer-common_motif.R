# find common motif between type 1,2,3
# homer
# ------------------------------------

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



### intersection
typeI_up <- union(homer_known[["typeI_up"]]$Motif.Name, homer_novo[["typeI_up"]]$Best.Match.Details)
typeII_up <- union(homer_known[["typeII_up"]]$Motif.Name, homer_novo[["typeII_up"]]$Best.Match.Details)
typeIII_up <- union(homer_known[["typeIII_up"]]$Motif.Name, homer_novo[["typeIII_up"]]$Best.Match.Details)

typeI_dn <- union(homer_known[["typeI_dn"]]$Motif.Name, homer_novo[["typeI_dn"]]$Best.Match.Details)
typeII_dn <- union(homer_known[["typeII_dn"]]$Motif.Name, homer_novo[["typeII_dn"]]$Best.Match.Details)
typeIII_dn <- union(homer_known[["typeIII_dn"]]$Motif.Name, homer_novo[["typeIII_dn"]]$Best.Match.Details)


common_up <- Reduce(intersect, list(typeI_up, typeII_up, typeIII_up))
common_dn <- Reduce(intersect, list(typeI_dn, typeII_dn, typeIII_dn))


## export
write.table(common_up, file = "homer-common_up.txt", quote = F, row.names = F, col.names = F)
write.table(common_dn, file = "homer-common_dn.txt", quote = F, row.names = F, col.names = F)