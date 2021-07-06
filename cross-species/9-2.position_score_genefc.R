# position score for pathways
# rank based on gene logFc
# ---------------------------

rm(list = ls())

library(tidyverse)


### DESeq2 results
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


## rank for genes
human_res$rank <- rank(human_res$log2FoldChange)
mouse_res$rank <- rank(mouse_res$log2FoldChange)


#### GO
#### -------------
library(topGO)


## topGO object
load("human_GOdata.Rdata")
load("mouse_GOdata.Rdata")


## topGO result
load("human_GOallres.Rdata")
load("mouse_GOallres.Rdata")


### position score
### s=2(R1 − R2)/n
### --------------


## postion score
# human

ps <- vector("list", 3)
names(ps) <- c("BP", "MF", "CC")

system.time({ for (go in names(human_allres)) {

        tmp <- c()

        selcGenes <- genesInTerm(human_GOdata[[go]], whichGO = human_allres[[go]]$GO.ID)

        for (term in names(selcGenes)) {

            de_inpath <- selcGenes[[term]]

            R1 <- human_res %>%
                    filter(id %in% de_inpath) %>%
                    dplyr::select(rank) %>% summarise_all(mean)

            R2 <- human_res %>%
                    filter(!id %in% de_inpath) %>%
                    dplyr::select(rank) %>% summarise_all(mean)

            S <- 2*(R1-R2) / nrow(human_res)

            tmp[term] <- S

            #ps[[go]] <- unlist(ps)

            #human_allres[[go]] <- cbind(human_allres[[go]], ps)
        }
        ps[[go]] <- unlist(tmp)
    }
})


human_bp <- cbind(human_allres[["BP"]], ps[["BP"]])
human_mf <- cbind(human_allres[["MF"]], ps[["MF"]])
human_cc <- cbind(human_allres[["CC"]], ps[["CC"]])

colnames(human_bp)[7] <- "ps"
colnames(human_mf)[7] <- "ps"
colnames(human_cc)[7] <- "ps"


# mouse

ps <- vector("list", 3)
names(ps) <- c("BP", "MF", "CC")

system.time({ for (go in names(mouse_allres)) {

        tmp <- c()

        selcGenes <- genesInTerm(mouse_GOdata[[go]], whichGO = mouse_allres[[go]]$GO.ID)

        for (term in names(selcGenes)) {

            de_inpath <- selcGenes[[term]]

            R1 <- mouse_res %>%
                    filter(id %in% de_inpath) %>%
                    dplyr::select(rank) %>% summarise_all(mean)

            R2 <- mouse_res %>%
                    filter(!id %in% de_inpath) %>%
                    dplyr::select(rank) %>% summarise_all(mean)

            S <- 2*(R1-R2) / nrow(mouse_res)

            tmp[term] <- S

        }
        ps[[go]] <- unlist(tmp)
    }
})


mouse_bp <- cbind(mouse_allres[["BP"]], ps[["BP"]])
mouse_mf <- cbind(mouse_allres[["MF"]], ps[["MF"]])
mouse_cc <- cbind(mouse_allres[["CC"]], ps[["CC"]])

colnames(mouse_bp)[7] <- "ps"
colnames(mouse_mf)[7] <- "ps"
colnames(mouse_cc)[7] <- "ps"


# save the GO res table with positionn score
res <- list(human_bp, human_mf, human_cc, mouse_bp, mouse_mf, mouse_cc)
names(res) <- c("human_bp", "human_mf", "human_cc", "mouse_bp", "mouse_mf", "mouse_cc")

for (go in names(res)) {
    write.table(res[[go]], file = paste0(go, "_gene_ps.txt"),
                quote = F, sep = "\t")
}


#### GSEA
#### -------------


## GSEA result
load("human_fgseaResTidy.Rdata")
load("mouse_fgseaResTidy.Rdata")



### position score
### s=2(R1 − R2)/n
### --------------


# human
ps <- c()

system.time({
    
    for (i in seq_len(nrow(h_fgseaResTidy))) {

    pathway <- h_fgseaResTidy$pathway[i]

    de_inpath <- h_fgseaResTidy$leadingEdge[[i]]

    R1 <- human_res %>%
                    filter(row %in% de_inpath) %>%
                    dplyr::select(rank) %>% summarise_all(mean)

    R2 <- human_res %>%
                    filter(!row %in% de_inpath) %>%
                    dplyr::select(rank) %>% summarise_all(mean)

    S <- 2*(R1-R2) / nrow(human_res)


    ps[i] <- S %>% unlist()
    names(ps)[i] <- pathway

    }
})

h_fgseaResTidy_gps <- cbind(h_fgseaResTidy, ps) %>% tibble()
save(h_fgseaResTidy_gps, file = "h_fgseaResTidy_gps.Rdata")


# mouse
ps <- c()

system.time({

    for (i in seq_len(nrow(m_fgseaResTidy))) {

    pathway <- m_fgseaResTidy$pathway[i]

    de_inpath <- m_fgseaResTidy$leadingEdge[[i]]

    R1 <- mouse_res %>%
                    filter(row %in% de_inpath) %>%
                    dplyr::select(rank) %>% summarise_all(mean)

    R2 <- mouse_res %>%
                    filter(!row %in% de_inpath) %>%
                    dplyr::select(rank) %>% summarise_all(mean)

    S <- 2*(R1-R2) / nrow(mouse_res)


    ps[i] <- S %>% unlist()
    names(ps)[i] <- pathway

    }
})

m_fgseaResTidy_gps <- cbind(m_fgseaResTidy, ps) %>% tibble()
save(m_fgseaResTidy_gps, file = "m_fgseaResTidy_gps.Rdata")
