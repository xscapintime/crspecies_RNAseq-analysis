# home de novo result
# dotplot
# ------------------------------------

rm(list = ls())

library(tidyverse)


# load table
load("denovo_motif_all.Rdata")

# change the name to human orthologs
dat_novo$ortholog <- toupper(dat_novo$name)

dat_novo$ortholog[grepl("BZIP", dat_novo$ortholog)] <- "bZIP_CREB"
dat_novo$ortholog[grepl("TATA", dat_novo$ortholog)] <- "TATA-Box"


# TF family, the solo one goes solo
regx <- c("^ARNT", "^ATF|^BATF", "^BCL", "^CEBP", "^DMRT", "^DUX", "^E2F", "^EOMES", "^ESR", "^ETS|^ETV", "^FOX", "^GATA", "^HIF",
            "^HNF", "^HOX", "^IRF", "KLF", "^MAF", "^MEF", "^MYC","MYB", "^NFAT", "^NKX", "^NR2", "^NRF", "^RAX", "^PRDM", "^RAR",
            "^RFX", "^RUNX", "^SIX", "^SMAD", "^SOX", "^SP", "^SPDEF", "^TBX", "^TCF", "^TEAD", "^ZBTB", "^ZFP", "^ZNF|^ZKSCAN|^ZSCAN")

fam <-  c("ARNT", "ATF", "BCL", "CEBP", "DMRT", "DUX", "E2F", "EOMES", "ESR", "ETS", "FOX", "GATA","HIF",
            "HNF", "HOX", "IRF", "KLF", "MAF", "MEF", "MYC", "MYB", "NFAT", "NKX", "NR2", "NRF", "RAX", "PRDM", "RAR",
            "RFX", "RUNX", "SIX", "SMAD", "SOX", "SP", "SPDEF", "TBX", "TCF", "TEAD", "ZBTB", "ZFP", "ZNF")

dat_novo$fam <- dat_novo$ortholog
for (i in 1:length(regx)) {
    dat_novo$fam[grepl(regx[i], dat_novo$ortholog)] <- fam[i]
}

# split 1 and 2,3
dat_novo$type <- unlist(lapply(X = as.character(dat_novo$group), FUN = function(x) {return(strsplit(x, split = "_")[[1]][1])}))
dat_novo$typegroup <- ifelse(dat_novo$type == "typeI", "type1", "type23")


# take the smallesr p-val one, if one group have duplicated motifs
dat_novo <- dat_novo %>% group_by(group, name) %>% filter(P.value == min(P.value))


# group by fam, oder by pval and percentage
tf_order <- (dat_novo %>% group_by(typegroup, fam) %>% arrange(P.value, desc(percent), .by_group = T))$ortholog %>% unique()
dat_novo$ortholog <- factor(dat_novo$ortholog, levels = tf_order)


# make TF family factors
fam_order <- (dat_novo %>% group_by(typegroup, fam) %>% summarise(mean(percent), min(P.value)) %>%
    arrange(`min(P.value)`, desc(`mean(percent)`), .by_group = T))$fam %>% unique()
dat_novo$fam <- factor(dat_novo$fam, levels = fam_order)





### plot ###
## colors
mycols <- c("#4B80DB", "#dc1a91","#dc1a84", "#dc1a64", "#DC1E1A")
morecols <- (grDevices::colorRampPalette(mycols))(125) ### expand the platte!


## all motif, for supp
theme_set(
    cowplot::theme_minimal_grid() +
    theme(axis.line  = element_line(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1 , size = 4),
        axis.text.y = element_text(size = 6),
        #axis.ticks = element_blank(),
        panel.border = element_blank(),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.key.size = unit(.3, "cm"),
        plot.margin = margin(t = 1.5, r = 1, b = .5, l = 1, unit = "cm")
    )
)


dat_novo %>% mutate(`P-value (-log10)` = -log10(P.value)) %>% #filter(percent > 10)

    ggplot(aes(x = ortholog, y = group, fill = `P-value (-log10)`, size = percent)) +
    geom_point(alpha = 1, shape = 21, colour = "black", stroke = .3) +
    ylab("") + xlab("") +
    scale_fill_gradientn(colours = morecols) +
    scale_y_discrete(limits = rev,
        labels = rev(c("Type I up", "Type II up", "Type III up", "Type I down", "Type II down", "Type III down"))) +
    labs(size = "% Motif in Target", fill = bquote(-log["10"] (P-value))) +
    cowplot::panel_border(color  = "black", size = .8)# +
    #scale_size(range = c(0, 3))


ggsave("motif_dotplot_novo_all.pdf", width = 16, height = 2.8)


## fig in the article
# FOX, ATF, MYC, MYB, RUNX, HNF, PPAR, CEBP, E2F, HIF and KLF

sel <- c("FOX", "ATF", "MYC", "MYB", "RUNX", "HNF", "PPAR", "CEBP", "E2F", "HIF", "KLF")

dat_sel <- dat_novo %>% filter(fam %in% sel)


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
        plot.margin = margin(t = 1.5, r = 1, b = .5, l = 1, unit = "cm")
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

ggsave("motif_dotplot_novo_sel.pdf", width = 7.5, height = 2.8)
