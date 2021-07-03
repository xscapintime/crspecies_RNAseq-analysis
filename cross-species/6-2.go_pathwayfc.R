# logFc for pathways
# ORA GO
# ------------------

rm(list = ls())
library(tidyverse)
library(topGO)


## Godata to retrive genes in GO term
load("human_GOdata.Rdata")
load("mouse_GOdata.Rdata")

load("human_GOallres.Rdata")
load("mouse_GOallres.Rdata")


## DE results
human_res <- read.table("human_deg.tsv")
mouse_res <- read.table("mouse_deg.tsv")

# back to ensembl id
index <- read.table("pairsidx.tsv")

human_res$id <- index$ensembl_gene_id.x[match(human_res$row, index$hgnc_symbol)]
mouse_res$id <- index$ensembl_gene_id.y[match(mouse_res$row, index$hgnc_symbol)]

# need remove NA pval/qval genes
human_res <- human_res %>% na.omit()
mouse_res <- mouse_res %>% na.omit()

# NA pval in mouse data, make them be in the same length
human_res <- human_res %>% filter(row %in% mouse_res$row)



##====== BP logFC ======##
## ----------------------
human_bp <- human_allres[["BP"]]
mouse_bp <- mouse_allres[["BP"]]


# check the "Annotated"
# turns out the "Annotated genes" all come from input genes
selcGenes <- genesInTerm(human_GOdata[["BP"]], whichGO = human_bp$GO.ID)

check <- c()
for (go in names(selcGenes)) {
    #print(go)
    tmp <- selcGenes[[go]] %in% human_res$id %>% table() %>% as.data.frame()
    if (nrow(tmp) != 1) {
        print(paste0(go, ": gene number doesn't fit"))
    } else {
            #print(paste0("this one good"))
            check[go] <- go
        }
    #check[go] <- selcGenes[[go]] %>% length()
}

length(check) == length(selcGenes)



## BP logFC
human_bp_fc <- list()

system.time({
    for (go in names(selcGenes)) {
        de_inpath <- selcGenes[[go]]
        
        fc <- human_res %>%
                filter(id %in% de_inpath) %>%
                dplyr::select(log2FoldChange) %>% summarise_all(mean)
        
        human_bp_fc[go] <- fc #%>% unlist()
}
})

human_bp_fc <- unlist(human_bp_fc)

write.table(human_bp_fc, file = "human_bp_fc.txt",
            quote = F, sep = "\t", row.names = T, col.names = F)



selcGenes <- genesInTerm(mouse_GOdata[["BP"]], whichGO = mouse_bp$GO.ID)

mouse_bp_fc <- list()

system.time({
    for (go in names(selcGenes)) {
        de_inpath <- selcGenes[[go]]
        
        fc <- mouse_res %>%
                filter(id %in% de_inpath) %>%
                dplyr::select(log2FoldChange) %>% summarise_all(mean)
        
        mouse_bp_fc[go] <- fc #%>% unlist()
}
})

mouse_bp_fc <- unlist(mouse_bp_fc)

write.table(mouse_bp_fc, file = "mouse_bp_fc.txt",
            quote = F, sep = "\t", row.names = T, col.names = F)



##====== MF logFC ======##
## -----------------------
human_mf <- human_allres[["MF"]]
mouse_mf <- mouse_allres[["MF"]]


# human
selcGenes <- genesInTerm(human_GOdata[["MF"]], whichGO = human_mf$GO.ID)

human_mf_fc <- list()

system.time({
    for (go in names(selcGenes)) {
        de_inpath <- selcGenes[[go]]
        
        fc <- human_res %>%
                filter(id %in% de_inpath) %>%
                dplyr::select(log2FoldChange) %>% summarise_all(mean)
        
        human_mf_fc[go] <- fc #%>% unlist()
}
})

human_mf_fc <- unlist(human_mf_fc)

write.table(human_mf_fc, file = "human_mf_fc.txt",
            quote = F, sep = "\t", row.names = T, col.names = F)


# mouse
selcGenes <- genesInTerm(mouse_GOdata[["MF"]], whichGO = mouse_mf$GO.ID)

mouse_mf_fc <- list()

system.time({
    for (go in names(selcGenes)) {
        de_inpath <- selcGenes[[go]]
        
        fc <- mouse_res %>%
                filter(id %in% de_inpath) %>%
                dplyr::select(log2FoldChange) %>% summarise_all(mean)
        
        mouse_mf_fc[go] <- fc #%>% unlist()
}
})

mouse_mf_fc <- unlist(mouse_mf_fc)

write.table(mouse_mf_fc, file = "mouse_mf_fc.txt",
            quote = F, sep = "\t", row.names = T, col.names = F)



##====== CC logFC ======##
## -----------------------
human_cc <- human_allres[["CC"]]
mouse_cc <- mouse_allres[["CC"]]


# human
selcGenes <- genesInTerm(human_GOdata[["CC"]], whichGO = human_cc$GO.ID)

human_cc_fc <- list()

system.time({
    for (go in names(selcGenes)) {
        de_inpath <- selcGenes[[go]]
        
        fc <- human_res %>%
                filter(id %in% de_inpath) %>%
                dplyr::select(log2FoldChange) %>% summarise_all(mean)
        
        human_cc_fc[go] <- fc #%>% unlist()
}
})

human_cc_fc <- unlist(human_cc_fc)

write.table(human_cc_fc, file = "human_cc_fc.txt",
            quote = F, sep = "\t", row.names = T, col.names = F)


# mouse
selcGenes <- genesInTerm(mouse_GOdata[["CC"]], whichGO = mouse_cc$GO.ID)

mouse_cc_fc <- list()

system.time({
    for (go in names(selcGenes)) {
        de_inpath <- selcGenes[[go]]
        
        fc <- mouse_res %>%
                filter(id %in% de_inpath) %>%
                dplyr::select(log2FoldChange) %>% summarise_all(mean)
        
        mouse_cc_fc[go] <- fc #%>% unlist()
}
})

mouse_cc_fc <- unlist(mouse_cc_fc)

write.table(mouse_cc_fc, file = "mouse_cc_fc.txt",
            quote = F, sep = "\t", row.names = T, col.names = F)