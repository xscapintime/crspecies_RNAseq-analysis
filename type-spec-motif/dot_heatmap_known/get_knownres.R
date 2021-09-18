# homer known result
# get motif name, p-val, percent
# need to merge up and dn
# ------------------------------------

rm(list = ls())

library(tidyverse)

## homer known
files_known <- list.files("../homer.v2", pattern = "knownResults.txt", recursive = T)

homer_known <- lapply(files_known, function(file) { read.csv(paste0("../homer.v2/", file), sep = "\t", header = T)})
homer_known <- lapply(homer_known, function(x) mutate(x, percent = as.numeric(sub("%","",`X..of.Target.Sequences.with.Motif`))) %>%
                        filter(P.value < 0.01) %>% select(1, 3, 10))

names(homer_known) <- unlist(lapply(X = files_known, FUN = function(x) {return(strsplit(x, split = "_MotifOutput")[[1]][1])}))

# group
homer_known <- Map(function(dat, grp) cbind(dat, group = grp), homer_known,
                unlist(lapply(X = names(homer_known), FUN = function(x) {return(x)})))

## reduce to table
dat_kn <- Reduce(rbind, homer_known)


## processing
# make group as factor
dat_kn$group <- factor(dat_kn$group, levels = c("typeI_up", "typeII_up", "typeIII_up", "typeI_dn", "typeII_dn", "typeIII_dn"))


# dealing with names
dat_kn$tf <- unlist(lapply(X = dat_kn$Motif.Name, FUN = function(x) {return(strsplit(x, split = "/")[[1]][1])}))
dat_kn$name <- unlist(lapply(X = dat_kn$tf, FUN = function(x) {return(strsplit(x, split = "\\(")[[1]][1])}))


# change the name to human orthologs
# some keep its original name
dat_kn$ortholog <- toupper(dat_kn$name)

dat_kn$ortholog[grepl("BHLHE40", dat_kn$ortholog)] <- "bHLHE40"
dat_kn$ortholog[grepl("BHLHE41", dat_kn$ortholog)] <- "bHLHE41"

dat_kn$ortholog[grepl("C-JUN-CRE", dat_kn$ortholog)] <- "c-JUN-CRE"

dat_kn$ortholog[grepl("C-MYC", dat_kn$ortholog)] <- "c-MYC"

dat_kn$ortholog[grepl("PIT1+1BP", dat_kn$ortholog)] <- "Pit1+1bp"

dat_kn$ortholog[grepl("RONIN", dat_kn$ortholog)] <- "Ronin"


# Tf family
regx <- c("MYB","^AP-1|^ATF|^BATF","^BACH","^bHLHE","^CDX","^CEBP","^CTCF","^DMRT","^E2F","^EGR","^ELF","^ELK","^ETS|^ETV","^EWS","^FOX","^FRA",
        "^GATA","^GFY","^HIF","^HOX","^HNF","^IRF","JUN","^KLF","^MEF","MYC","^NFY","^NFKB","^NKX","^NPAS","^NRF","^PBX","^PHOX2A",
        "^PIT1","^PITX1","^PRDM","^PPAR","^RAR","^RBPJ","^RFX","^SIX","^SMAD","^SOX","^SP\\d","^SPI|^PU.1","^STAT","^TEAD","^USF","^ZBTB","^ZFP","^ZNF|^ZSCAN|^ZKSCAN")

fam <- c("MYB","AP-1/ATF","BACH","bHLHE","CDX","CEBP","CTCF","DMRT","E2F","EGR","ELF","ELK","ETS","EWS","FOX","FRA","GATA","GFY","HIF","HOX","HNF",
        "IRF","JUN","KLF","MEF","MYC","NFY","NFKB","NKX","NPAS","NRF","PBX","PHOX2A","PIT1","PITX1","PRDM","PPAR","RAR","RBPJ","RFX","SIX","SMAD","SOX",
        "SP","SPI","STAT","TEAD","USF","ZBTB","ZFP","ZNF")

dat_kn$fam <- NA
for (i in 1:length(regx)) {
    dat_kn$fam[grepl(regx[i], dat_kn$ortholog)] <- fam[i]
}

dat_kn$fam <- ifelse(is.na(dat_kn$fam), dat_kn$ortholog, dat_kn$fam)


# split 1 and 2,3
dat_kn$type <- unlist(lapply(X = as.character(dat_kn$group), FUN = function(x) {return(strsplit(x, split = "_")[[1]][1])}))
dat_kn$type <- factor(dat_kn$type, levels = c("typeI", "typeII", "typeIII"))

dat_kn$typegroup <- ifelse(dat_kn$type == "typeI", "type1", "type23")


# merge type up dn to type deg
dat_kn <- dat_kn %>% group_by(type, ortholog) %>% filter(P.value == min(P.value))

# take the smallest p-val one, if one group have duplicated motifs
# dat_kn <- dat_kn %>% group_by(group, name) %>% filter(P.value == min(P.value))


## order
# family order, by mean of percentage and the smallest p-val in the family
# make TF family factors
fam_order <- (dat_kn %>% group_by(typegroup, fam) %>% summarise(mean(percent), min(P.value)) %>%
    arrange(`min(P.value)`, desc(`mean(percent)`), .by_group = T))$fam %>% unique()
dat_kn$fam <- factor(dat_kn$fam, levels = fam_order)

# TF order within a family
# group by typegroup and fam, order by pval and percentage
tf_order <- (dat_kn %>% group_by(fam) %>%
    arrange(P.value, mean(percent), .by_group = T))$ortholog %>% unique()
dat_kn$ortholog <- factor(dat_kn$ortholog, levels = tf_order)



### save the table
save(dat_kn, file = "known_motif_all.Rdata")
write.table(dat_kn, file = "known_motif_all.txt", quote = F, sep = "\t", row.names = F, col.names = T)
