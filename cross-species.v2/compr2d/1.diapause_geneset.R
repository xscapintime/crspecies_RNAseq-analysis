# get mouse diapause gene set from DESeq2 reult
# compare with the one in Duy et al., 2021
# they also used Boroviak et al., 2015 data
# -----------------------------------------------


rm(list = rm())

library(tidyverse)

## mouse deg table
mouse_res_t2 <- read.table("../diff-expr/mouse_deg_t2.tsv")
mouse_res_t2 <- mouse_res_t2 %>% na.omit()

mouse_res_t3 <- read.table("../diff-expr/mouse_deg_t3.tsv")
mouse_res_t3 <- mouse_res_t3 %>% na.omit()


## sig genes
# abs(logFC) > 2 and adj p-val < 0.05

# up
diapuse_up_t2 <- mouse_res_t2 %>% filter(log2FoldChange > 2 & padj < 0.05) %>%
                    dplyr::select(row)

diapuse_up_t3 <- mouse_res_t3 %>% filter(log2FoldChange > 2 & padj < 0.05) %>%
                    dplyr::select(row)


# down
diapuse_down_t2 <- mouse_res_t2 %>% filter(log2FoldChange < -2 & padj < 0.05) %>%
                    dplyr::select(row)

diapuse_down_t3 <- mouse_res_t3 %>% filter(log2FoldChange < -2 & padj < 0.05) %>%
                    dplyr::select(row)


## compare with gene set in Duy 2021 cancer discovery paper
# just load it
library(readxl)
duy <- read_excel("../../cross-species/diapause_set/252204_2_supp_6841527_qmhm08.xlsx")

duy_up <- duy$DIAPAUSE_UP_BOROVIAK %>% na.omit
duy_up <- duy_up[-1]

duy_down <- duy$BOROVIAK_DIAPAUSE_DN %>% na.omit
duy_down <- duy_down[-1]

# up to up, down to down
diapuse_up_t3$row %in% duy_up %>% table
diapuse_down_t3$row %in% duy_down %>% table

mouse_res_t2 %>% filter(row %in% duy_up) %>% dim
mouse_res_t2 %>% filter(row %in% duy_down) %>% dim


# turn to gmt-like file
diapuse_set_t2 <- list(diapuse_up_t2$row, diapuse_down_t2$row, duy_up, duy_down)
diapuse_set_t3 <- list(diapuse_up_t3$row, diapuse_down_t3$row, duy_up, duy_down)

names(diapuse_set_t2) <- c("Diapause_up_t2", "Diapause_dn_t2", "Diapause_up_Duy2021", "Diapause_dn_Duy2021")
names(diapuse_set_t3) <- c("Diapause_up_t3", "Diapause_dn_t3", "Diapause_up_Duy2021", "Diapause_dn_Duy2021")

save(diapuse_set_t2, file = "diapuse_set_t2.Rdata")
save(diapuse_set_t3, file = "diapuse_set_t3.Rdata")