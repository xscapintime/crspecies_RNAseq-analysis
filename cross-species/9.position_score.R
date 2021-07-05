# position score for pathways
# rank based on pathway logFc
# ---------------------------

rm(list = ls())


library(tidyverse)


#### GO
#### -------------


## pahtway logFC
human_bp_fc <- read.table("human_bp_fc.txt", row.names = 1)
human_mf_fc <- read.table("human_mf_fc.txt", row.names = 1)
human_cc_fc <- read.table("human_cc_fc.txt", row.names = 1)

mouse_bp_fc <- read.table("mouse_bp_fc.txt", row.names = 1)
mouse_mf_fc <- read.table("mouse_mf_fc.txt", row.names = 1)
mouse_cc_fc <- read.table("mouse_cc_fc.txt", row.names = 1)


## topGO result
load("human_GOallres.Rdata")
load("mouse_GOallres.Rdata")


# BP
human_bp <- cbind(human_allres[["BP"]], human_bp_fc) %>% arrange(desc(V2))
human_mf <- cbind(human_allres[["MF"]], human_mf_fc) %>% arrange(desc(V2))
human_cc <- cbind(human_allres[["CC"]], human_cc_fc) %>% arrange(desc(V2))

mouse_bp <- cbind(mouse_allres[["BP"]], mouse_bp_fc) %>% arrange(desc(V2))
mouse_mf <- cbind(mouse_allres[["MF"]], mouse_mf_fc) %>% arrange(desc(V2))
mouse_cc <- cbind(mouse_allres[["CC"]], mouse_cc_fc) %>% arrange(desc(V2))



### position score
### s=2(R1 − R2)/n
### --------------


# rank
human_bp$rank <- rank(human_bp$V2)
human_mf$rank <- rank(human_mf$V2)
human_cc$rank <- rank(human_cc$V2)

mouse_bp$rank <- rank(mouse_bp$V2)
mouse_mf$rank <- rank(mouse_mf$V2)
mouse_cc$rank <- rank(mouse_cc$V2)


# postion score
res <- list(human_bp, human_mf, human_cc, mouse_bp, mouse_mf, mouse_cc)
names(res) <- c("human_bp", "human_mf", "human_cc", "mouse_bp", "mouse_mf", "mouse_cc")

for (go in names(res)) {

    tmp <- res[[go]]

    for (j in 1:nrow(tmp)) {
        R1 <- tmp$rank[j]
        R2 <- mean(tmp$rank[-j])
        tmp$ps[j] <- 2*(R1-R2)/nrow(tmp)
        
        res[[go]]$ps <- tmp$ps
    }
}


# human_bp <- res$human_bp
# human_mf <- res$human_mf
# human_cc <- res$human_cc

# mouse_bp <- res$mouse_bp
# mouse_mf <- res$mouse_mf
# mouse_cc <- res$mouse_cc


# save the GO res table with positionn score
for (go in names(res)) {
    write.table(res[[go]], file = paste0(go, "_fc_ps.txt"),
                quote = F, sep = "\t")
}



#### GSEA
#### -------------