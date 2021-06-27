#library(org.Mm.eg.db)


rm(list = ls())
options(scipen = 999)
library(tidyverse)

#BiocManager::install("topGO")
library(topGO)
library(biomaRt)


## deg table
human_res <- read.table("human_deg.tsv")
mouse_res <- read.table("mouse_deg.tsv")

## ensembl id and symble
human_symbol <- read.table("human_idsyb.tsv")
mouse_symbol_tohuman <- read.table("mouse_idsyb_mapped2hugo.tsv")

# get mouse symbol
mouse <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
mouse_symbol <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol"),
                    filters = "ensembl_gene_id",
                    values = mouse_symbol_tohuman$ensembl_gene_id, mart = mouse)
mouse_symbol_df <- inner_join(mouse_symbol, mouse_symbol_tohuman,
                                by = c("ensembl_gene_id" = "ensembl_gene_id"))
mouse_symbol_df  <- mouse_symbol_df %>% na_if("") %>% na.omit() %>% distinct()


## join deg table and symbol
human_res_syb <- inner_join(human_res, human_symbol,
                            by = c("row" = "ensembl_gene_id"))

mouse_res_syb <- inner_join(mouse_res, mouse_symbol_df,
                            by = c("row" = "ensembl_gene_id"))


## gene lists
# universe lists
human_res_syb$de <- ifelse(abs(human_res_syb$log2FoldChange) >= 2 &
                            human_res_syb$padj < 0.01, TRUE, FALSE)
human_universe <- as.numeric(human_res_syb$de) %>% factor() %>%
        setNames(human_res_syb$hgnc_symbol)

mouse_res_syb$de <- ifelse(abs(mouse_res_syb$log2FoldChange) >= 0.8 &
                            mouse_res_syb$padj < 0.05, TRUE, FALSE)
mouse_universe <- as.numeric(mouse_res_syb$de) %>% factor() %>%
        setNames(mouse_res_syb$mgi_symbol)


## GO annotation
# set ensembl datasets
human <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
#mouse <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")

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