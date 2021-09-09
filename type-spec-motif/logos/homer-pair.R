# for homer result, only one up and dn motif are 3-types0-common
# find pairwise-overlapped motifs
# ---------------------------------------------------------------

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
# typeI_up_idx <- union(homer_known[["typeI_up"]]$Motif.Name, homer_novo[["typeI_up"]]$Best.Match.Details)
# typeII_up_idx <- union(homer_known[["typeII_up"]]$Motif.Name, homer_novo[["typeII_up"]]$Best.Match.Details)
# typeIII_up_idx <- union(homer_known[["typeIII_up"]]$Motif.Name, homer_novo[["typeIII_up"]]$Best.Match.Details)

# typeI_dn_idx <- union(homer_known[["typeI_dn"]]$Motif.Name, homer_novo[["typeI_dn"]]$Best.Match.Details)
# typeII_dn_idx <- union(homer_known[["typeII_dn"]]$Motif.Name, homer_novo[["typeII_dn"]]$Best.Match.Details)
# typeIII_dn_idx <- union(homer_known[["typeIII_dn"]]$Motif.Name, homer_novo[["typeIII_dn"]]$Best.Match.Details)


typeI_up <- full_join(homer_known[["typeI_up"]][c("Motif.Name","P.value")],
                    homer_novo[["typeI_up"]][c("Best.Match.Details", "P.value")],
                    by = c("Motif.Name" = "Best.Match.Details"), suffix = c(".known", ".novo"))

typeII_up <- full_join(homer_known[["typeII_up"]][c("Motif.Name","P.value")],
                    homer_novo[["typeII_up"]][c("Best.Match.Details", "P.value")],
                    by = c("Motif.Name" = "Best.Match.Details"), suffix = c(".known", ".novo"))

typeIII_up <- full_join(homer_known[["typeIII_up"]][c("Motif.Name","P.value")],
                    homer_novo[["typeIII_up"]][c("Best.Match.Details", "P.value")],
                    by = c("Motif.Name" = "Best.Match.Details"), suffix = c(".known", ".novo"))



typeI_dn <- full_join(homer_known[["typeI_dn"]][c("Motif.Name","P.value")],
                    homer_novo[["typeI_dn"]][c("Best.Match.Details", "P.value")],
                    by = c("Motif.Name" = "Best.Match.Details"), suffix = c(".known", ".novo"))

typeII_dn <- full_join(homer_known[["typeII_dn"]][c("Motif.Name","P.value")],
                    homer_novo[["typeII_dn"]][c("Best.Match.Details", "P.value")],
                    by = c("Motif.Name" = "Best.Match.Details"), suffix = c(".known", ".novo"))

typeIII_dn <- full_join(homer_known[["typeIII_dn"]][c("Motif.Name","P.value")],
                    homer_novo[["typeIII_dn"]][c("Best.Match.Details", "P.value")],
                    by = c("Motif.Name" = "Best.Match.Details"), suffix = c(".known", ".novo"))



### type I vs II
typeIvsII_up <- inner_join(typeI_up, typeII_up, by = c("Motif.Name" = "Motif.Name"), suffix = c(".I", ".II"))
typeIvsII_dn <- inner_join(typeI_dn, typeII_dn, by = c("Motif.Name" = "Motif.Name"), suffix = c(".I", ".II"))

### type I vs III
typeIvsIII_up <- inner_join(typeI_up, typeIII_up, by = c("Motif.Name" = "Motif.Name"), suffix = c(".I", ".III"))
typeIvsIII_dn <- inner_join(typeI_dn, typeIII_dn, by = c("Motif.Name" = "Motif.Name"), suffix = c(".I", ".III"))

### type II vs III
typeIIvsIII_up <- inner_join(typeII_up, typeIII_up, by = c("Motif.Name" = "Motif.Name"), suffix = c(".II", ".III"))
typeIIvsIII_dn <- inner_join(typeII_dn, typeIII_dn, by = c("Motif.Name" = "Motif.Name"), suffix = c(".II", ".III"))


### to exclude the 1-2-3-common motif
all_pair <- list(typeIvsII_up, typeIvsII_dn, typeIvsIII_up, typeIvsIII_dn,
                typeIIvsIII_up, typeIIvsIII_dn)
names(all_pair) <- c("typeIvsII_up", "typeIvsII_dn", "typeIvsIII_up", "typeIvsIII_dn",
                    "typeIIvsIII_up", "typeIIvsIII_dn")

# common up, Nur77(NR)/K562-NR4A1-ChIP-Seq(GSE31363)/Homer
all_pair[c(1,3,5)] <- lapply(all_pair[c(1,3,5)] , function(x) filter(x, Motif.Name != "Nur77(NR)/K562-NR4A1-ChIP-Seq(GSE31363)/Homer"))

# common dn, ZKSCAN1(Zf)/HepG2-ZKSCAN1-ChIP-Seq(Encode)/Homer
all_pair[c(2,4,6)] <- lapply(all_pair[c(2,4,6)] , function(x) filter(x, Motif.Name != "ZKSCAN1(Zf)/HepG2-ZKSCAN1-ChIP-Seq(Encode)/Homer"))


### export
for (i in 1:length(all_pair)) {
    write.table(all_pair[[i]], file = paste0("homer-", names(all_pair[i]),"_fin.txt"), quote = F, row.names = F, col.names = F, sep = "\t")
}