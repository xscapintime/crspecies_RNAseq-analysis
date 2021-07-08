# calculating postion score based on gene rank
# s=2(R1 − R2)/n
# ---------------------------------------------


rm(list = ls())

library(tidyverse)


## gene ranks
load("ranks.Rdata")


## genes mapped to gene sets
load("gene2sets_10cut.Rdata")



### position score based on gene rank
## for each pathway

# human axis position
ps <- c()

system.time({

    for (i in 1:length(map_fin)) {

    pathway <- names(map_fin)[i]

    de_inpath <- map_fin[[i]] %>% as.logical
    names(de_inpath) <- names(map_fin[[i]])

    R1 <- human_rank[de_inpath] %>% mean()

    R2 <- human_rank[!de_inpath] %>% mean()

    S <- 2*(R1-R2) / length(human_rank)  # just for gene number


    ps[i] <- S
    names(ps)[i] <- pathway

    }
})

human_ps <- ps


# mouse axis position
ps <- c()

system.time({

    for (i in 1:length(map_fin)) {

    pathway <- names(map_fin)[i]

    de_inpath <- map_fin[[i]] %>% as.logical
    names(de_inpath) <- names(map_fin[[i]])

    R1 <- mouse_rank[de_inpath] %>% mean()

    R2 <- mouse_rank[!de_inpath] %>% mean()

    S <- 2*(R1-R2) / length(mouse_rank)  # just for gene number


    ps[i] <- S
    names(ps)[i] <- pathway

    }
})

mouse_ps <- ps


save(human_ps, mouse_ps, file = "position_score.Rdata")