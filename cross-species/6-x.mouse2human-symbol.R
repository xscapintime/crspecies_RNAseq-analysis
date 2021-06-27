#####==== mouse ensembl id to huamn id to GO id ====#####
# feel wrong
# ----------


mouse_2_human_go_id <- getLDS(attributes = c("ensembl_gene_id"),
                            filters = "ensembl_gene_id",
                            values =  mouse_res$row, mart = mouse,
                            attributesL = c("ensembl_gene_id", "go_id"), martL = human)


mouse_2_human_symb <- getLDS(attributes = c("ensembl_gene_id"),
                            filters = "ensembl_gene_id",
                            values =  mouse_res$row, mart = mouse,
                            attributesL = c("ensembl_gene_id"), martL = human)

mouse_2_human_go_id <- getBM(attributes = c("ensembl_gene_id", "go_id", "namespace_1003"),
                            filters = "ensembl_gene_id",
                            values = mouse_2_human_symb$Gene.stable.ID.1, mart = human)

mouse_2_human_go_id <- inner_join(mouse_2_human_go_id, mouse_2_human_symb,
                                by = c("ensembl_gene_id" = "Gene.stable.ID.1"))

tmp <- inner_join(mouse_2_human_go_id, mouse_go_id,
                by = c("Gene.stable.ID" = "ensembl_gene_id"))

unique(tmp$go_id.x) %>% length()
unique(tmp$go_id.y) %>% length()



#### if entrez id
mouse_go_id3 <- getBM(attributes = c("entrezgene_id","ensembl_gene_id", "go_id", "namespace_1003"),
                    filters = "ensembl_gene_id",
                    values = mouse_res$row, mart = mouse)
mouse_go_id3 <- mouse_go_id3 %>% na_if("") %>% na.omit() %>% distinct()

unique(mouse_go_id3$ensembl_gene_id) %>% length
unique(mouse_go_id3$entrezgene_id) %>% length








### no use?
## ensembl id and symble
human_symbol <- read.table("human_idsyb.tsv")
mouse_symbol_tohuman <- read.table("mouse_idsyb_mapped2hugo.tsv")


## join deg table and symbol
human_res_syb <- inner_join(human_res, human_symbol,
                            by = c("row" = "ensembl_gene_id"))

mouse_res_syb <- inner_join(mouse_res, mouse_symbol_tohuman,
                            by = c("row" = "ensembl_gene_id"))
