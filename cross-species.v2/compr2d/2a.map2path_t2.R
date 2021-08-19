# map all genes to GO, KEGG and diapause set
# use gmtPathways fucntion from fgsea package
# for type 2 data
# -------------------------------------------


rm(list = ls())

library(tidyverse)

### gene sets
## load GO and KEGG gmt
library(fgsea)

bp <- gmtPathways("../../cross-species/gmtdata/c5.go.bp.v7.4.symbols.gmt")
cc <- gmtPathways("../../cross-species/gmtdata/c5.go.cc.v7.4.symbols.gmt")
mf <- gmtPathways("../../cross-species/gmtdata/c5.go.mf.v7.4.symbols.gmt")

kegg <- gmtPathways("../../cross-species/gmtdata/c2.cp.kegg.v7.4.symbols.gmt")


## diapause gene set
load("diapuse_set_t2.Rdata")

diapuse_set <- diapuse_set_t2 # easy name


### DEG table
## DESeq2 results
human_res <- read.table("../diff-expr/human_deg_t2.tsv")
mouse_res <- read.table("../diff-expr/mouse_deg_t2.tsv")

# need remove NA pval/qval genes
human_res <- human_res %>% na.omit()
mouse_res <- mouse_res %>% na.omit()

# NA pval in mouse data, make them be in the same length
idx <- intersect(human_res$row, mouse_res$row)

human_res <- human_res %>% filter(row %in% idx)
mouse_res <- mouse_res %>% filter(row %in% idx)

setdiff(human_res$row, mouse_res$row)


## rank based on logFC
human_rank <- setNames(rank(human_res$log2FoldChange), human_res$row)
mouse_rank <- setNames(rank(mouse_res$log2FoldChange), mouse_res$row)

save(human_rank, mouse_rank, file = "ranks_t2.Rdata")


### map all genes to the gene sets

sets <- list(bp, mf, cc, kegg, diapuse_set)
names(sets) <- c("bp", "mf", "cc", "kegg", "diapuse_set")


map <- list()

system.time({

    for (set in names(sets)) {

        mapset <- list()
        cat <- sets[[set]]

        for (i in 1:length(cat)) {
            tmp <- names(human_rank) %in% cat[[i]] %>% as.numeric()
            names(tmp) <- names(human_rank)   # just for gene names
            mapset[[i]] <- tmp
            names(mapset)[i] <- names(cat)[i]
        }

    map[[set]] <- mapset

    }

})

save(map, file = "gene2sets_t2.Rdata")


## remove pathway <10 input genes
# calculate the sum

thesum <- list()

system.time({
    for (set in names(map)) {
        cat <- map[[set]]

        tmp <- c()
        for (i in 1:length(cat)) {
            tmp[i] <- sum(cat[[i]])
            names(tmp)[i] <- names(cat)[i]
        }

    thesum[[set]] <- tmp

    }
})

thecheck <- c(unlist(thesum$bp), unlist(thesum$mf), unlist(thesum$cc),
                unlist(thesum$kegg), unlist(thesum$diapuse_set))

torm <- thecheck[thecheck < 10]


# install.packages("rlist")
library(rlist)

map_fin <- c(list.remove(map$bp, names(torm)),
            list.remove(map$mf, names(torm)),
            list.remove(map$cc, names(torm)),
            list.remove(map$kegg, names(torm)),
            list.remove(map$diapuse_set, names(torm)))


save(map_fin, file = "gene2sets_10cut_t2.Rdata")