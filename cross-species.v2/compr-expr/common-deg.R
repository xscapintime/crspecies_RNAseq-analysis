# find common DEG between human and mouse data
# t2 and t3
# FC > 2 and p-val < 0.05
# ---------------------------------------------------


rm(list = ls())

library(tidyverse)


## type 2
human_res_t2 <- read.table("../diff-expr/human_deg_t2.tsv")
mouse_res_t2 <- read.table("../diff-expr/mouse_deg_t2.tsv")

# up
human_up_t2 <- human_res_t2 %>% filter(log2FoldChange > 2 & padj < 0.05) %>% dplyr::select(row)
mouse_up_t2 <- mouse_res_t2 %>% filter(log2FoldChange > 2 & padj < 0.05) %>% dplyr::select(row)

com_up_t2 <- intersect(human_up_t2, mouse_up_t2)

# down
human_dn_t2 <- human_res_t2 %>% filter(log2FoldChange < -2 & padj < 0.05) %>% dplyr::select(row)
mouse_dn_t2 <- mouse_res_t2 %>% filter(log2FoldChange < -2 & padj < 0.05) %>% dplyr::select(row)

com_dn_t2 <- intersect(human_dn_t2, mouse_dn_t2)


## type 3
human_res_t3 <- read.table("../diff-expr/human_deg_t3.tsv")
mouse_res_t3 <- read.table("../diff-expr/mouse_deg_t3.tsv")

# up
human_up_t3 <- human_res_t3 %>% filter(log2FoldChange > 2 & padj < 0.05) %>% dplyr::select(row)
mouse_up_t3 <- mouse_res_t3 %>% filter(log2FoldChange > 2 & padj < 0.05) %>% dplyr::select(row)

com_up_t3 <- intersect(human_up_t3, mouse_up_t3)

# down
human_dn_t3 <- human_res_t3 %>% filter(log2FoldChange < -2 & padj < 0.05) %>% dplyr::select(row)
mouse_dn_t3 <- mouse_res_t3 %>% filter(log2FoldChange < -2 & padj < 0.05) %>% dplyr::select(row)

com_dn_t3 <- intersect(human_dn_t3, mouse_dn_t3)



## export
write.table(com_up_t2, file = "common_up_deg_t2.tsv", quote = F, sep = "\t", row.names = F, col.names = F)
write.table(com_dn_t2, file = "common_dn_deg_t2.tsv", quote = F, sep = "\t", row.names = F, col.names = F)

write.table(com_up_t3, file = "common_up_deg_t3.tsv", quote = F, sep = "\t", row.names = F, col.names = F)
write.table(com_dn_t3, file = "common_dn_deg_t3.tsv", quote = F, sep = "\t", row.names = F, col.names = F)