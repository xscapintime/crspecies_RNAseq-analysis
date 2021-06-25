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



mouse_glist <- mouse_res$padj
names(mouse_glist) <- mouse_res$row

mouse_GOdata_2 <- new("topGOdata",
                        ontology = "BP",
                        allGenes = mouse_glist,
                        annot = annFUN.org,
                        geneSel = topDiffGenes,
                        nodeSize = 10,
                        mapping = "org.Mm.eg",
                        ID = "ensembl")








