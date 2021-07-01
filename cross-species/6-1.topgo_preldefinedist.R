# Predefined list of interesting gene
# -------------ORA v1 ---------------

rm(list = ls())
#options(scipen = 999)
library(tidyverse)

#BiocManager::install("topGO")
library(topGO)

# annotated by biomart
# --------------------
library(biomaRt)


## deg table
human_res <- read.table("human_deg.tsv")
mouse_res <- read.table("mouse_deg.tsv")


# back to ensembl id
index <- read.table("pairsidx.tsv")

human_res$id <- index$ensembl_gene_id.x[match(human_res$row, index$hgnc_symbol)]
mouse_res$id <- index$ensembl_gene_id.y[match(mouse_res$row, index$hgnc_symbol)]


# need remove NA pval/qval genes
human_res <- human_res %>% na.omit()
mouse_res <- mouse_res %>% na.omit() ## mouse data has NA

# NA pval in mouse data, make them be in the same length
human_res <- human_res %>% filter(row %in% mouse_res$row)


## gene lists
# universe lists
human_res$de <- ifelse(human_res$padj < 0.01, TRUE, FALSE)
human_universe <- as.numeric(human_res$de) %>% factor() %>%
        setNames(human_res$id)

mouse_res$de <- ifelse(mouse_res$padj < 0.01, TRUE, FALSE)
mouse_universe <- as.numeric(mouse_res$de) %>% factor() %>%
        setNames(mouse_res$id)


## GO annotation
# set ensembl datasets
human <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
mouse <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")

# GO id
human_go_id <- getBM(attributes = c("ensembl_gene_id", "go_id", "namespace_1003"),
                        filters = "ensembl_gene_id",
                        values = human_res$id, mart = human)
human_go_id <- human_go_id %>% na_if("") %>% na.omit() %>% distinct()

mouse_go_id <- getBM(attributes = c("ensembl_gene_id", "go_id", "namespace_1003"),
                     filters = "ensembl_gene_id",
                     values = mouse_res$id, mart = mouse)
mouse_go_id <- mouse_go_id %>% na_if("") %>% na.omit() %>% distinct()


### topGO
### -----
## build gene2GO list
h_gene2GO <- unstack(human_go_id[, c(2, 1)])

m_gene2GO <- unstack(mouse_go_id[, c(2, 1)])


## topgo object and fisher test
# BP

go <- c("BP", "MF", "CC")

human_GOdata <- list()
human_allres <- list()

for (p in go) {


        human_GOdata[[p]] <- new("topGOdata",
                                ontology = p,
                                description = "human, biomart, predefined genelist",
                                allGenes = human_universe,
                                annot = annFUN.gene2GO,
                                gene2GO = h_gene2GO,
                                nodeSize = 10)

        human_res_fisher <- runTest(human_GOdata[[p]], algorithm = "classic", statistic = "fisher")

        human_allres[[p]] <- GenTable(human_GOdata[[p]], classicFisher = human_res_fisher,
                                        orderBy = "classicFisher",
                                        ranksOf = "classicFisher",
                                        topNodes = length(usedGO(human_GOdata[[p]])))

}

save(human_GOdata, file = "human_GOdata.Rdata")
save(human_allres, file = "human_GOallres.Rdata")

write.table(human_allres[["BP"]], file = "human_predfres_bp.tsv", quote = F, sep = "\t")
write.table(human_allres[["MF"]], file = "human_predfres_mf.tsv", quote = F, sep = "\t")
write.table(human_allres[["CC"]], file = "human_predfres_cc.tsv", quote = F, sep = "\t")



mouse_GOdata <- list()
mouse_allres <- list()

for (p in go) {

        mouse_GOdata[[p]] <- new("topGOdata",
                                ontology = p,
                                description = "mouse, biomart, predefined genelist",
                                allGenes = mouse_universe,
                                annot = annFUN.gene2GO,
                                gene2GO = m_gene2GO,
                                nodeSize = 10)

        mouse_res_fisher <- runTest(mouse_GOdata[[p]], algorithm = "classic", statistic = "fisher")

        mouse_allres[[p]] <- GenTable(mouse_GOdata[[p]], classicFisher = mouse_res_fisher,
                                        orderBy = "classicFisher",
                                        ranksOf = "classicFisher",
                                        topNodes = length(usedGO(mouse_GOdata[[p]])))
}

save(mouse_GOdata, file = "mouse_GOdata.Rdata")
save(mouse_allres, file = "mouse_GOallres.Rdata")

write.table(mouse_allres[["BP"]], file = "mouse_predfres_bp.tsv", quote = F, sep = "\t")
write.table(mouse_allres[["MF"]], file = "mouse_predfres_mf.tsv", quote = F, sep = "\t")
write.table(mouse_allres[["CC"]], file = "mouse_predfres_cc.tsv", quote = F, sep = "\t")
