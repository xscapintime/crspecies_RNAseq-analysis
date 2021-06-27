# Predefined list of interesting gene
# -------------ORA v1 ---------------

rm(list = ls())
options(scipen = 999)
library(tidyverse)

#BiocManager::install("topGO")
library(topGO)

# annotated by biomart
# --------------------
library(biomaRt)


## deg table
human_res <- read.table("human_deg.tsv")
mouse_res <- read.table("mouse_deg.tsv")

# need remove NA pval/qval genes
human_res <- human_res %>% na.omit()
mouse_res <- mouse_res %>% na.omit() ## mouse data has NA

## gene lists
# universe lists
human_res$de <- ifelse(abs(human_res$log2FoldChange) >= 2 &
                        human_res$padj < 0.01, TRUE, FALSE)
human_universe <- as.numeric(human_res$de) %>% factor() %>%
        setNames(human_res$row)

mouse_res$de <- ifelse(abs(mouse_res$log2FoldChange) >= 2 &
                        mouse_res$padj < 0.01, TRUE, FALSE)
mouse_universe <- as.numeric(mouse_res$de) %>% factor() %>%
        setNames(mouse_res$row)


## GO annotation
# set ensembl datasets
human <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
mouse <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")

# GO id
human_go_id <- getBM(attributes = c("ensembl_gene_id", "go_id", "namespace_1003"),
                        filters = "ensembl_gene_id",
                        values = human_res$row, mart = human)
human_go_id <- human_go_id %>% na_if("") %>% na.omit() %>% distinct()

mouse_go_id <- getBM(attributes = c("ensembl_gene_id", "go_id", "namespace_1003"),
                     filters = "ensembl_gene_id",
                     values = mouse_res$row, mart = mouse)
mouse_go_id <- mouse_go_id %>% na_if("") %>% na.omit() %>% distinct()


### topGO
### -----
## build gene2GO list
h_gene2GO <- unstack(human_go_id[, c(2, 1)])

m_gene2GO <- unstack(mouse_go_id[, c(2, 1)])


## topgo object and fisher test
# BP
human_GOdata <- new("topGOdata",
                    ontology = "BP",
                    description = "human, biomart, predefined genelist",
                    allGenes = human_universe,
                    annot = annFUN.gene2GO,
                    gene2GO = h_gene2GO)

mouse_GOdata <- new("topGOdata",
                    ontology = "BP",
                    description = "mouse, biomart, predefined genelist",
                    allGenes = mouse_universe,
                    annot = annFUN.gene2GO,
                    gene2GO = m_gene2GO)


human_res_fisher <- runTest(human_GOdata, algorithm = "classic", statistic = "fisher")
human_allres_bp <- GenTable(human_GOdata, classicFisher = human_res_fisher,
                            orderBy = "classicFisher",
                            ranksOf = "classicFisher",
                            topNodes = length(usedGO(human_GOdata)))

mouse_res_fisher <- runTest(mouse_GOdata, algorithm = "classic", statistic = "fisher")
mouse_allres_bp <- GenTable(mouse_GOdata, classicFisher = mouse_res_fisher,
                            orderBy = "classicFisher",
                            ranksOf = "classicFisher",
                            topNodes = length(usedGO(mouse_GOdata)))

write.table(human_allres_bp, file = "human_predfres_bp.tsv", quote = F, sep = "\t")
write.table(mouse_allres_bp, file = "mouse_predfres_bp.tsv", quote = F, sep = "\t")


# MF
human_GOdata <- new("topGOdata",
                    ontology = "MF",
                    description = "human, biomart, predefined genelist",
                    allGenes = human_universe,
                    annot = annFUN.gene2GO,
                    gene2GO = h_gene2GO)

mouse_GOdata <- new("topGOdata",
                    ontology = "MF",
                    description = "mouse, biomart, predefined genelist",
                    allGenes = mouse_universe,
                    annot = annFUN.gene2GO,
                    gene2GO = m_gene2GO)


human_res_fisher <- runTest(human_GOdata, algorithm = "classic", statistic = "fisher")
human_allres_mf <- GenTable(human_GOdata, classicFisher = human_res_fisher,
                            orderBy = "classicFisher",
                            ranksOf = "classicFisher",
                            topNodes = length(usedGO(human_GOdata)))

mouse_res_fisher <- runTest(mouse_GOdata, algorithm = "classic", statistic = "fisher")
mouse_allres_mf <- GenTable(mouse_GOdata, classicFisher = mouse_res_fisher,
                            orderBy = "classicFisher",
                            ranksOf = "classicFisher",
                            topNodes = length(usedGO(mouse_GOdata)))

write.table(human_allres_mf, file = "human_predfres_mf.tsv", quote = F, sep = "\t")
write.table(mouse_allres_mf, file = "mouse_predfres_mf.tsv", quote = F, sep = "\t")



# CC
human_GOdata <- new("topGOdata",
                    ontology = "CC",
                    description = "human, biomart, predefined genelist",
                    allGenes = human_universe,
                    annot = annFUN.gene2GO,
                    gene2GO = h_gene2GO)

mouse_GOdata <- new("topGOdata",
                    ontology = "CC",
                    description = "mouse, biomart, predefined genelist",
                    allGenes = mouse_universe,
                    annot = annFUN.gene2GO,
                    gene2GO = m_gene2GO)


human_res_fisher <- runTest(human_GOdata, algorithm = "classic", statistic = "fisher")
human_allres_cc <- GenTable(human_GOdata, classicFisher = human_res_fisher,
                            orderBy = "classicFisher",
                            ranksOf = "classicFisher",
                            topNodes = length(usedGO(human_GOdata)))

mouse_res_fisher <- runTest(mouse_GOdata, algorithm = "classic", statistic = "fisher")
mouse_allres_cc <- GenTable(mouse_GOdata, classicFisher = mouse_res_fisher,
                            orderBy = "classicFisher",
                            ranksOf = "classicFisher",
                            topNodes = length(usedGO(mouse_GOdata)))

write.table(human_allres_cc, file = "human_predfres_cc.tsv", quote = F, sep = "\t")
write.table(mouse_allres_cc, file = "mouse_predfres_cc.tsv", quote = F, sep = "\t")