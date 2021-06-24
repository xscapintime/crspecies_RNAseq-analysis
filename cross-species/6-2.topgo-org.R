rm(list = ls())
options(scipen = 999)
library(tidyverse)

#BiocManager::install("topGO")
library(topGO)
library(org.Hs.eg.db)
library(org.Mm.eg.db)


## deg table
human_res <- read.table("human_deg.tsv")
mouse_res <- read.table("mouse_deg.tsv")


## gene lists
# universe lists
human_res$de <- ifelse(abs(human_res$log2FoldChange) >= 2 &
                            human_res$padj < 0.01, TRUE, FALSE)
human_universe <- as.numeric(human_res$de) %>% factor() %>%
        setNames(human_res$row)

mouse_res$de <- ifelse(mouse_res$padj < 0.05, TRUE, FALSE)
mouse_universe <- as.numeric(mouse_res$de) %>% factor() %>%
        setNames(mouse_res$row)


## topgo object and fisher test
# BP
human_GOdata <- new("topGOdata",
                    ontology = "BP",
                    allGenes = human_universe,
                    annot = annFUN.org,
                    mapping = "org.Hs.eg",
                    ID = "ensembl")

mouse_GOdata <- new("topGOdata",
                    ontology = "BP",
                    allGenes = mouse_universe,
                    annot = annFUN.org,
                    mapping = "org.Mm.eg",
                    ID = "ensembl")





mouse_res_fisher <- runTest(mouse_GOdata, algorithm = "classic", statistic = "fisher")
mouse_allres_bp <- GenTable(mouse_GOdata, classicFisher = mouse_res_fisher,
                            orderBy = "classicFisher",
                            ranksOf = "classicFisher",
                            topNodes = length(usedGO(mouse_GOdata)))





### if using gene score?

topDiffGenes <- function(allScore) {
        return(allScore < 0.01)
        }



mouse_glist <- human_res$padj
names(mouse_glist) <- mouse_res$row

mouse_GOdata_2 <- new("topGOdata",
                        ontology = "BP",
                        allGenes = mouse_glist,
                        annot = annFUN.org,
                        geneSel = topDiffGenes,
                        nodeSize = 10,
                        mapping = "org.Mm.eg",
                        ID = "ensembl")








## GO annotation
# set ensembl datasets
human <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
mouse <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")

# GO id
human_go_id <- getBM(attributes = c("hgnc_symbol", "go_id", "namespace_1003"),
                    filters = "hgnc_symbol",
                    values = human_res_syb$hgnc_symbol, mart = human)
human_go_id <- human_go_id %>% na_if("") %>% na.omit() %>% distinct()

mouse_go_id <- getBM(attributes = c("mgi_symbol", "go_id", "namespace_1003"),
                    filters = "mgi_symbol",
                    values = mouse_res_syb$mgi_symbol, mart = mouse)
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
                    allGenes = human_universe,
                    annot = annFUN.gene2GO,
                    gene2GO = h_gene2GO)

mouse_GOdata <- new("topGOdata",
                    ontology = "BP",
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