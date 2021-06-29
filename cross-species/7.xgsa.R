# XGSA
# ----

## looks like need to compare human DE with mouse unverse GO
## or can I use topGO abtained GO list objects?


rm(list = ls())
library(tidyverse)

#devtools::install_github('VCCRI/XGSA')
#options(download.file.method = "wininet")
#options(download.file.method = "libcurl")

# or maybe solutions in https://github.com/r-lib/remotes/issues/130

library(xgsa)

## huamn DE
# huamn DE table
human_res <- read.table("human_deg.tsv")

# select some genes
human_de_genes <- human_res %>% filter(padj < 0.01 & log2FoldChange > 2) %>% dplyr::select(row) %>% deframe()

# gene universe, ensembl id and symbol
human_ensembl_symbol_map <- get_ENSEMBL_symbol_map(species = "hsapiens")


## XGSA object
human_xgsa <- new_XGSA_dataset(species = "hsapiens",
                                 data = list(humanDegenes = human_de_genes),
                                 type = 'genesetlist', name = 'humanDegenes',
                                 universe = unique(human_ensembl_symbol_map$ensembl_gene_id))


## mouse GO (all) to against with
mouse_bp <- get_GO("mmusculus", ontologies = "biological_process")
mouse_bp <- mouse_bp[lapply(mouse_bp, length) > 10 & lapply(mouse_bp, length) < 500]
mouse_bp_xgsa <- new_XGSA_dataset(species = "mmusculus",
                                    data = mouse_bp, type = 'genesetlist',
                                    name = "mouseBP", universe = unique(unlist(mouse_bp)))


## run XGSA!
human_vs_mousebp <- run_XGSA_test(human_xgsa, mouse_bp_xgsa)

resulting.pvals <- lapply(human_vs_mousebp, function(X){ X[["pvals"]] })
resulting.overlap.genes <- lapply(human_vs_mousebp, function(X){ X[["genes"]] })

adjusted.pvals <- p.adjust(unlist(resulting.pvals), method = "BH")
names(adjusted.pvals) <- unlist(lapply(strsplit(names(adjusted.pvals) ,"\\."), function(X){return(X[[2]])}))

mouse_bp_names <- get_GO_names("mmusculus")
names(adjusted.pvals) <- mouse_bp_names[match(names(adjusted.pvals), mouse_bp_names$go_id), "name_1006"]


## finally the result
significant.GO.Terms <- adjusted.pvals[which(adjusted.pvals < 0.05)]
