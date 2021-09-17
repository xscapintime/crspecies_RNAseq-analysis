# find common motif between type 1,2,3
# put homer known and de novo output seprately?
# ------------------------------------

rm(list = ls())

library(tidyverse)


### homer res
## list files
files_known <- list.files("../homer.v2", pattern = "knownResults.txt", recursive = T)
files_novo <- list.files("../homer.v2", pattern = "homerResults.txt", recursive = T)


### homer result
## homer known
homer_known <- lapply(files_known, function(file) { read.csv(paste0("../homer.v2/", file), sep = "\t", header = T) %>% filter(P.value < 0.01)})
homer_known <- lapply(homer_known, function(x) mutate(x, percent = as.numeric(sub("%","",`X..of.Target.Sequences.with.Motif`))))
homer_known <- lapply(homer_known, function(x) filter(x, P.value < 0.01))# & percent >= 10))
names(homer_known) <- unlist(lapply(X = files_known, FUN = function(x) {return(strsplit(x, split = "_MotifOutput")[[1]][1])}))


## homer de novo
homer_novo <- lapply(files_novo, function(file) { read.csv(paste0("../homer.v2/", file), sep = "\t", header = T) %>% filter(P.value < 0.01)})
homer_novo <- lapply(homer_novo, function(x) mutate(x, percent = as.numeric(sub("%","",`X..of.Targets`))))
homer_novo <- lapply(homer_novo, function(x) filter(x, P.value < 0.01)) #& score > 0.8))# & percent >= 10))
names(homer_novo) <- unlist(lapply(X = files_novo, FUN = function(x) {return(strsplit(x, split = "_MotifOutput")[[1]][1])}))


### intersection
typeI_up <- full_join(homer_known[["typeI_up"]][c("Motif.Name","P.value","percent")],
                    homer_novo[["typeI_up"]][c("Best.Match.Details", "P.value","percent", "score")],
                    by = c("Motif.Name" = "Best.Match.Details"), suffix = c(".known", ".novo"))

typeII_up <- full_join(homer_known[["typeII_up"]][c("Motif.Name","P.value","percent")],
                    homer_novo[["typeII_up"]][c("Best.Match.Details", "P.value","percent", "score")],
                    by = c("Motif.Name" = "Best.Match.Details"), suffix = c(".known", ".novo"))

typeIII_up <- full_join(homer_known[["typeIII_up"]][c("Motif.Name","P.value","percent")],
                    homer_novo[["typeIII_up"]][c("Best.Match.Details", "P.value","percent", "score")],
                    by = c("Motif.Name" = "Best.Match.Details"), suffix = c(".known", ".novo"))


typeI_dn <- full_join(homer_known[["typeI_dn"]][c("Motif.Name","P.value","percent")],
                    homer_novo[["typeI_dn"]][c("Best.Match.Details", "P.value","percent", "score")],
                    by = c("Motif.Name" = "Best.Match.Details"), suffix = c(".known", ".novo"))

typeII_dn <- full_join(homer_known[["typeII_dn"]][c("Motif.Name","P.value","percent")],
                    homer_novo[["typeII_dn"]][c("Best.Match.Details", "P.value","percent", "score")],
                    by = c("Motif.Name" = "Best.Match.Details"), suffix = c(".known", ".novo"))

typeIII_dn <- full_join(homer_known[["typeIII_dn"]][c("Motif.Name","P.value","percent")],
                    homer_novo[["typeIII_dn"]][c("Best.Match.Details", "P.value","percent", "score")],
                    by = c("Motif.Name" = "Best.Match.Details"), suffix = c(".known", ".novo"))



common_up <- Reduce(intersect, list(typeI_up$Motif.Name, typeII_up$Motif.Name, typeIII_up$Motif.Name))
common_dn <- Reduce(intersect, list(typeI_dn$Motif.Name, typeII_dn$Motif.Name, typeIII_dn$Motif.Name))


## export
write.table(common_up, file = "homer-common_up.txt", quote = F, row.names = F, col.names = F)
write.table(common_dn, file = "homer-common_dn.txt", quote = F, row.names = F, col.names = F)


allspec <- list(typeI_up, typeII_up, typeIII_up, typeI_dn, typeII_dn, typeIII_dn)
names(allspec) <- c("typeI_up", "typeII_up", "typeIII_up", "typeI_dn", "typeII_dn", "typeIII_dn")

for (i in 1:length(allspec)) {
    write.table(allspec[[i]], file = paste0("homer-", names(allspec[i]),".txt"), quote = F, row.names = F, col.names = T, sep = "\t")
}
