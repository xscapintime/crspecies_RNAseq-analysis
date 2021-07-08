# get mouse diapause gene set from DESeq2 reult
# compare with the one in Duy et al., 2021
# they also used Boroviak et al., 2015 data
# -----------------------------------------------


rm(list = rm())

library(tidyverse)

## mouse deg table
mouse_res <- read.table("../mouse_deg.tsv")

mouse_res <- mouse_res %>% na.omit()


## sig genes
# abs(logFC) > 2 and adj p-val < 0.05

# up
diapuse_up <- mouse_res %>% filter(log2FoldChange > 2 & padj < 0.05) %>%
                    dplyr::select(row)



# down
diapuse_down <- mouse_res %>% filter(log2FoldChange < -2 & padj < 0.05) %>%
                    dplyr::select(row)



## compare with gene set in Duy 2021 cancer discovery paper
library(readxl)
duy <- read_excel("../diapause_set/252204_2_supp_6841527_qmhm08.xlsx")

duy_up <- duy$DIAPAUSE_UP_BOROVIAK %>% na.omit
duy_up <- duy_up[-1]

duy_down <- duy$BOROVIAK_DIAPAUSE_DN %>% na.omit
duy_down <- duy_down[-1]

# up to up, down to down
diapuse_up$row %in% duy_up %>% table
diapuse_down$row %in% duy_down %>% table

mouse_res %>% filter(row %in% duy_up) %>% dim
mouse_res %>% filter(row %in% duy_down) %>% dim


# turn to gmt-like file
diapuse_set <- list(diapuse_up$row, diapuse_down$row, duy_up, duy_down)
names(diapuse_set) <- c("Diapause_up", "Diapause_dn", "Diapause_up_Duy2021", "Diapause_dn_Duy2021")

save(diapuse_set, file = "diapuse_set.Rdata")