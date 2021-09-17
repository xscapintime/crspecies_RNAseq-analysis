# homer de novo result
# get motif name, p-val, percent, and score
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


## processing
# make group as factor
dat_novo$group <- factor(dat_novo$group, levels = c("typeI_up", "typeII_up", "typeIII_up", "typeI_dn", "typeII_dn", "typeIII_dn"))


# split 1 and 2,3
dat_novo$type <- unlist(lapply(X = as.character(dat_novo$group), FUN = function(x) {return(strsplit(x, split = "_")[[1]][1])}))
dat_novo$typegroup <- ifelse(dat_novo$type == "typeI", "type1", "type23")


# dealing with names
dat_novo$tf <- unlist(lapply(X = dat_novo$Best.Match.Details, FUN = function(x) {return(strsplit(x, split = "/")[[1]][1])}))
dat_novo$name <- unlist(lapply(X = dat_novo$tf, FUN = function(x) {return(strsplit(x, split = "\\(")[[1]][1])}))
dat_novo$name[grepl("_", dat_novo$name)] <- unlist(lapply(X = dat_novo$name[grepl("_", dat_novo$name)], FUN = function(x) {return(strsplit(x, split = "1_")[[1]][2])}))


# change the name to human orthologs
dat_novo$ortholog <- toupper(dat_novo$name)

dat_novo$ortholog[grepl("BZIP", dat_novo$ortholog)] <- "bZIP_CREB"
dat_novo$ortholog[grepl("TATA", dat_novo$ortholog)] <- "TATA-Box"


# TF family, the solo one goes solo
regx <- c("^ARNT", "^ATF|^BATF", "^BCL", "^CEBP", "^DMRT", "^DUX", "^E2F", "^EOMES", "^ESR", "^ETS|^ETV", "^FOX", "^GATA", "^HIF",
            "^HNF", "^HOX", "^IRF", "KLF", "^MAF", "^MEF", "^MYC","MYB", "^NFAT", "^NKX", "^NR2", "^NRF", "^RAX", "^PRDM", "^RAR",
            "^RFX", "^RUNX", "^SIX", "^SMAD", "^SOX", "^SP", "^SPDEF", "^TBX", "^TCF", "^TEAD", "^ZBTB", "^ZFP", "^ZNF|^ZKSCAN|^ZSCAN")

fam <-  c("ARNT", "ATF", "BCL", "CEBP", "DMRT", "DUX", "E2F", "EOMES", "ESR", "ETS", "FOX", "GATA","HIF",
            "HNF", "HOX", "IRF", "KLF", "MAF", "MEF", "MYC", "MYB", "NFAT", "NKX", "NR2", "NRF", "RAX", "PRDM", "RAR",
            "RFX", "RUNX", "SIX", "SMAD", "SOX", "SP", "SPDEF", "TBX", "TCF", "TEAD", "ZBTB", "ZFP", "ZNF")

dat_novo$fam <- dat_novo$ortholog
for (i in 1:length(regx)) {
    dat_novo$fam[grepl(regx[i], dat_novo$ortholog)] <- fam[i]
}


# take the smallest p-val one, if one group have duplicated motifs
dat_novo <- dat_novo %>% group_by(group, name) %>% filter(P.value == min(P.value))


## order
# TF order within a family
# group by typegroup and fam, order by pval and percentage
tf_order <- (dat_novo %>% group_by(typegroup, fam) %>% arrange(P.value, desc(percent), .by_group = T))$ortholog %>% unique()
dat_novo$ortholog <- factor(dat_novo$ortholog, levels = tf_order)


# family order, by mean of percentage and the smallest p-val in the family
# make TF family factors
fam_order <- (dat_novo %>% group_by(typegroup, fam) %>% summarise(mean(percent), min(P.value)) %>%
    arrange(`min(P.value)`, desc(`mean(percent)`), .by_group = T))$fam %>% unique()
dat_novo$fam <- factor(dat_novo$fam, levels = fam_order)



### save the table
save(dat_novo, file = "denovo_motif_all.Rdata")
write.table(dat_novo, file = "denovo_motif_all.txt", quote = F, sep = "\t", row.names = F, col.names = T)