# find common motif between type 1,2,3
# and specific ones
# meme
# ------------------------------------

rm(list = ls())

library(tidyverse)


### streme-tomtom res
## list files
files <- list.files("../streme/tomtom", pattern = "tomtom.tsv", recursive = T)


### tomtom result
tomtom <- lapply(files, function(file) { read.csv(paste0("../streme/tomtom/", file), sep = "\t", header = T) %>% filter(q.value < 0.01)})
names(tomtom) <- unlist(lapply(X = files, FUN = function(x) {return(strsplit(x, split = "_tomtom")[[1]][1])}))


### common
## maybe need q.value of each
common_up <- Reduce(intersect, list(tomtom[["typeI_up"]]$Target_ID, tomtom[["typeII_up"]]$Target_ID, tomtom[["typeIII_up"]]$Target_ID))
common_dn <- Reduce(intersect, list(tomtom[["typeI_dn"]]$Target_ID, tomtom[["typeII_dn"]]$Target_ID, tomtom[["typeIII_dn"]]$Target_ID))


## export
write.table(common_up, file = "meme-common_up.txt", quote = F, row.names = F, col.names = F)
write.table(common_dn, file = "meme-common_dn.txt", quote = F, row.names = F, col.names = F)


### specific motif
typeI_spec_up <-  tomtom[["typeI_up"]] %>% filter(Target_ID %in% setdiff(tomtom[["typeI_up"]]$Target_ID, common_up))
typeII_spec_up <- setdiff(tomtom[["typeII_up"]]$Target_ID, common_up)
typeIII_spec_up <- setdiff(tomtom[["typeIII_up"]]$Target_ID, common_up)

typeI_spec_dn <- setdiff(tomtom[["typeI_dn"]]$Target_ID, common_dn)
typeII_spec_dn <- setdiff(tomtom[["typeII_dn"]]$Target_ID, common_dn)
typeIII_spec_dn <- setdiff(tomtom[["typeIII_dn"]]$Target_ID, common_dn)


allspec <- list(typeI_spec_up, typeII_spec_up, typeIII_spec_up, typeI_spec_dn, typeII_spec_dn, typeIII_spec_dn)
names(allspec) <- c("typeI_spec_up", "typeII_spec_up", "typeIII_spec_up", "typeI_spec_dn", "typeII_spec_dn", "typeIII_spec_dn")

for (i in 1:length(allspec)) {
    write.table(allspec[[i]], file = paste0("meme-", names(allspec[i]),".txt"), quote = F, row.names = F, col.names = F)
}
