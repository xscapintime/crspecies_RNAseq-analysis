# home de novo result
# get motif name, p-val, percet, and score
# ------------------------------------

rm(list = ls())

library(tidyverse)

## homer de novo
files_novo <- list.files("../homer.v2", pattern = "homerResults.txt", recursive = T)

homer_novo <- lapply(files_novo, function(file) { read.csv(paste0("../homer.v2/", file), sep = "\t", header = T)})
homer_novo <- lapply(homer_novo, function(x) mutate(x, percent = as.numeric(sub("%","",`X..of.Targets`))))
homer_novo <- lapply(homer_novo, function(x) filter(x, P.value < 0.01) %>% select("Best.Match.Details", "P.value", "percent", "score")) #& score > 0.8))# & percent >= 10))
names(homer_novo) <- unlist(lapply(X = files_novo, FUN = function(x) {return(strsplit(x, split = "_MotifOutput")[[1]][1])}))

# group
homer_novo <- Map(function(dat, grp) cbind(dat, group = grp), homer_novo,
                unlist(lapply(X = names(homer_novo), FUN = function(x) {return(x)})))

## reduce to table
dat_novo <- Reduce(rbind, homer_novo)

# make group as factor
dat_novo$group <- factor(dat_novo$group, levels = c("typeI_up", "typeII_up", "typeIII_up", "typeI_dn", "typeII_dn", "typeIII_dn"))

# dealing with ames
dat_novo$tf <- unlist(lapply(X = dat_novo$Best.Match.Details, FUN = function(x) {return(strsplit(x, split = "/")[[1]][1])}))
dat_novo$name <- unlist(lapply(X = dat_novo$tf, FUN = function(x) {return(strsplit(x, split = "\\(")[[1]][1])}))
dat_novo$name[grepl("_", dat_novo$name)] <- unlist(lapply(X = dat_novo$name[grepl("_", dat_novo$name)], FUN = function(x) {return(strsplit(x, split = "1_")[[1]][2])}))


## save the table
save(dat_novo, file = "denovo_motif_all.Rdata")
write.table(dat_novo, file = "denovo_motif_all.txt", quote = F, sep = "\t", row.names = F, col.names = T)