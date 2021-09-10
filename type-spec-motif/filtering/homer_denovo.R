# find common motif between type 1,2,3
# homer de novo
# intersect by 
# ------------------------------------

rm(list = ls())

library(tidyverse)

files_novo <- list.files("../homer", pattern = "homerResults.txt", recursive = T)

## homer de novo
homer_novo <- lapply(files_novo, function(file) { read.csv(paste0("../homer/", file), sep = "\t", header = T) %>% filter(P.value < 0.01)})
names(homer_novo) <- unlist(lapply(X = files_novo, FUN = function(x) {return(strsplit(x, split = "_MotifOutput")[[1]][1])}))


## up
common_up <- Reduce(intersect, list(homer_novo[["typeI_up"]]$Motif, homer_novo[["typeII_up"]]$Motif, homer_novo[["typeIII_up"]]$Motif))
common_dn <- Reduce(intersect, list(homer_novo[["typeI_dn"]]$Motif, homer_novo[["typeII_dn"]]$Motif, homer_novo[["typeIII_dn"]]$Motif))

#========== no overlap =============#