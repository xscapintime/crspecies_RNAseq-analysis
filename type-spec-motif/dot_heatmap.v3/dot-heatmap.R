# dot heatmap
# looks like @malikPluripotencyReprogrammingCompetent2019 Fig 2b
# ---------------------------------------------------------------

rm(list = ls())

library(tidyverse)

### load homer tables
load("motif_all.Rdata")


#### ==================================================== ####
#### fig in the article
# FOX, ATF, MYC, MYB, RUNX, HNF, PPAR, CEBP, E2F, HIF and KLF

dat_sel <- dat %>% filter(across(Motif.Name, ~ grepl("^FOX|^Fox|^ATF|^Atf|^MYC|^Myc|^MYB|^Myb|^RUNX|^Runx|^HNF|^Hnf|^PPAR|^Ppar|^CEBP|^Cebp|^E2F|^E2f|^HIF|^Hif|^KLF|^Flf", .)))
dat_sel$name <- unlist(lapply(X = dat_sel$tf, FUN = function(x) {return(strsplit(x, split = "\\(")[[1]][1])}))

dat_sel$name[grepl("\\(", dat_sel$tf)] <- unlist(lapply(X = dat_sel$tf[grepl("\\(", dat_sel$tf)], FUN = function(x) {return(strsplit(x, split = "\\(")[[1]][1])}))
dat_sel$name[!grepl("\\(", dat_sel$tf)] <- dat_sel$tf[!grepl("\\(", dat_sel$tf)]

# if same TF and same group, take the one with higher percentage
dat_sel <- dat_sel %>% group_by(group, name) %>% filter(percent == max(percent))

# Tf family
sel <- c("^FOX|^Fox", "^ATF|^Atf", "^MYC|^Myc","^MYB|^Myb","^RUNX|^Runx","^HNF|^Hnf","^PPAR|^Ppar","^CEBP|^Cebp","^E2F|^E2f","^HIF|^Hif","^KLF|^Flf")
fam <- c("FOX",  "ATF", "MYC", "MYB", "RUNX", "HNF", "PPAR", "CEBP", "E2F", "HIF", "KLF")

dat_sel$fam <- NA
for (i in 1:length(sel)) {
    dat_sel$fam[grepl(sel[i], dat_sel$name)] <- fam[i]
}

# de
dat_sel$de <- unlist(lapply(X = as.character(dat_sel$group), FUN = function(x) {return(strsplit(x, split = "_")[[1]][2])}))

# make TF family factors
percent_mean <- (dat_sel %>% group_by(fam) %>% summarise(mean(percent), min(P.value)) %>% arrange(`min(P.value)`, desc(`mean(percent)`)))$fam
dat_sel$fam <- factor(dat_sel$fam, levels = percent_mean)

# change the name to human orthologs
dat_sel$ortholog <- toupper(dat_sel$name)

# oder by pval and percentage
tf_order <- (dat_sel %>% group_by(fam, group) %>% arrange(desc(percent), P.value, .by_group = T))$ortholog %>% unique()
dat_sel$ortholog <- factor(dat_sel$ortholog, levels = tf_order)



## plot ==================================================================================
## colors
mycols <- c("#4B80DB", "#dc1a91","#dc1a84", "#dc1a64", "#DC1E1A")
morecols <- (grDevices::colorRampPalette(mycols))(125) ### expand the platte!


# selected
theme_set(
    cowplot::theme_minimal_grid() +
    theme(axis.line  = element_line(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1 , size = 8),
        axis.text.y = element_text(size = 8),
        #axis.ticks = element_blank(),
        panel.border = element_blank(),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.key.size = unit(.3, "cm"),
        plot.margin = margin(t = 1, r = 1, b = .5, l = 1, unit = "cm")
    )
)

dat_sel %>% mutate(`P-value (-log10)` = -log10(P.value)) %>% #filter(percent > 10 & `P-value (-log10)` > 3) %>%
    ggplot(aes(x = ortholog, y = group, fill = `P-value (-log10)`, size = percent)) +
    geom_point(alpha = 1, shape = 21, colour = "black", stroke = .5) +
    ylab("") + xlab("") +
    scale_fill_gradientn(colours = morecols) +
    #scale_size(trans = "reverse") +
    scale_y_discrete(limits = rev,
        labels = rev(c("Type I up", "Type II up", "Type III up", "Type I down", "Type II down", "Type III down"))) +
    labs(size = "% Motif in Target", fill = bquote(-log["10"] (P-value))) +
    cowplot::panel_border(color  = "black", size = .8)

ggsave("motif_dotplot_sel.pdf", width = 10.5, height = 2.8)