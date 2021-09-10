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

## add q.value
common_up_df <- cbind(tomtom[["typeI_up"]] %>% filter(Target_ID %in% common_up) %>% group_by(Target_ID) %>% summarise(min(q.value)) %>% select(`min(q.value)`),
                        tomtom[["typeII_up"]] %>% filter(Target_ID %in% common_up) %>% group_by(Target_ID) %>% summarise(min(q.value)) %>% select(`min(q.value)`),
                        tomtom[["typeIII_up"]] %>% filter(Target_ID %in% common_up) %>% group_by(Target_ID) %>% summarise(min(q.value)) %>% select(`min(q.value)`))
dimnames(common_up_df) <- list(sort(common_up), names(tomtom)[c(2,4,6)])
common_up_df <- common_up_df[order(common_up_df),] %>% na.omit



common_dn_df <- cbind(tomtom[["typeI_dn"]] %>% filter(Target_ID %in% common_dn) %>% group_by(Target_ID) %>% summarise(min(q.value)) %>% select(`min(q.value)`),
                        tomtom[["typeII_dn"]] %>% filter(Target_ID %in% common_dn) %>% group_by(Target_ID) %>% summarise(min(q.value)) %>% select(`min(q.value)`),
                        tomtom[["typeIII_dn"]] %>% filter(Target_ID %in% common_dn) %>% group_by(Target_ID) %>% summarise(min(q.value)) %>% select(`min(q.value)`))
dimnames(common_dn_df) <- list(sort(common_dn), names(tomtom)[c(1,3,5)])
common_dn_df <- common_dn_df[order(common_dn_df),] %>% na.omit


## export
write.table(common_up_df, file = "meme-common_up.txt", quote = F, row.names = T, col.names = T)
write.table(common_dn_df, file = "meme-common_dn.txt", quote = F, row.names = T, col.names = T)


### specific motif
typeI_spec_up <-  tomtom[["typeI_up"]] %>% filter(Target_ID %in% setdiff(tomtom[["typeI_up"]]$Target_ID, common_up)) %>% group_by(Target_ID) %>% summarise(min(q.value))
typeII_spec_up <- tomtom[["typeII_up"]] %>% filter(Target_ID %in% setdiff(tomtom[["typeII_up"]]$Target_ID, common_up)) %>% group_by(Target_ID) %>% summarise(min(q.value))
typeIII_spec_up <- tomtom[["typeIII_up"]] %>% filter(Target_ID %in% setdiff(tomtom[["typeIII_up"]]$Target_ID, common_up)) %>% group_by(Target_ID) %>% summarise(min(q.value))

# sort by q-val
typeI_spec_up <- typeI_spec_up[order(typeI_spec_up$`min(q.value)`),] %>% na.omit
typeII_spec_up <- typeII_spec_up[order(typeII_spec_up$`min(q.value)`),] %>% na.omit
typeIII_spec_up <- typeIII_spec_up[order(typeIII_spec_up$`min(q.value)`),] %>% na.omit


typeI_spec_dn <- tomtom[["typeI_dn"]] %>% filter(Target_ID %in% setdiff(tomtom[["typeI_dn"]]$Target_ID, common_dn)) %>% group_by(Target_ID) %>% summarise(min(q.value))
typeII_spec_dn <- tomtom[["typeII_dn"]] %>% filter(Target_ID %in% setdiff(tomtom[["typeII_dn"]]$Target_ID, common_dn)) %>% group_by(Target_ID) %>% summarise(min(q.value))
typeIII_spec_dn <- tomtom[["typeIII_dn"]] %>% filter(Target_ID %in% setdiff(tomtom[["typeIII_dn"]]$Target_ID, common_dn)) %>% group_by(Target_ID) %>% summarise(min(q.value))

# sort by q-val
typeI_spec_dn <- typeI_spec_dn[order(typeI_spec_dn$`min(q.value)`),] %>% na.omit
typeII_spec_dn <- typeII_spec_dn[order(typeII_spec_dn$`min(q.value)`),] %>% na.omit
typeIII_spec_dn <- typeIII_spec_dn[order(typeIII_spec_dn$`min(q.value)`),] %>% na.omit


allspec <- list(typeI_spec_up, typeII_spec_up, typeIII_spec_up, typeI_spec_dn, typeII_spec_dn, typeIII_spec_dn)
names(allspec) <- c("typeI_spec_up", "typeII_spec_up", "typeIII_spec_up", "typeI_spec_dn", "typeII_spec_dn", "typeIII_spec_dn")

for (i in 1:length(allspec)) {
    write.table(allspec[[i]], file = paste0("meme-", names(allspec[i]),".txt"), quote = F, row.names = F, col.names = F)
}
