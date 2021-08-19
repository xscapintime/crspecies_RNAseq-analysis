# rank-MANOVA for each gene set
# for type 2 data
# -----------------------------


rm(list = ls())

library(tidyverse)


## gene ranks
load("ranks_t2.Rdata")


## genes mapped to gene sets
load("gene2sets_10cut_t2.Rdata")


## MANOVA Hotelling’s T2 test
# make an empty data frame

# number of all gene sets
# set_len <- c()
# for (set in names(map)) {
#     set_len[set] <- length(map[[set]])
# }
# set_len <- sum(unlist(set_len))


dat <- cbind(human_rank, mouse_rank) %>% data.frame

pval_all <- c()

system.time({

    for (i in 1:length(map_fin)) {
        dat$group <- map_fin[[i]] %>% factor()
        # MANOVA
        fit <- manova(cbind(human_rank, mouse_rank) ~ group, data = dat)
        summ <- summary(fit, test = "Hotelling-Lawley")
        pval  <- summ$stats[11]
        pval_all[i] <- pval
        }

})

names(pval_all) <- names(map_fin)

save(pval_all, file = "manova_pval_t2.Rdata")

# another way to get p-val from manova result
# https://stackoverflow.com/a/64074660/14498100
#summary(m, "Hotelling-Lawley")$stats["Species", "Pr(>F)"]



## ajust p-val
bh_adjpval <- p.adjust(pval_all, method = "BH")
save(bh_adjpval, file = "manova_bh_adjpval_t2.Rdata")
