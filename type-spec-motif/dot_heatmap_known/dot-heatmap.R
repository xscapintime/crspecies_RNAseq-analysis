# dot heatmap
# looks like @malikPluripotencyReprogrammingCompetent2019 Fig 2b
# ---------------------------------------------------------------

rm(list = ls())

library(tidyverse)

### load homer tables
load("motif_known_all.Rdata")


#### ==================================================== ####
#### fig in the article
# FOX, ATF, MYC, MYB, RUNX, HNF, PPAR, CEBP, E2F, HIF and KLF
# sel <- c("FOX", "^AP-1", "MYC", "MYB", "RUNX", "HNF4a", "PPAR", "CEBP", "E2F", "HIF", "KLF")

dat_sel <- dat %>% filter(across(fam, ~ grepl("^FOX|ATF|^MYC|MYB|^RUNX|^HNF|^PPAR|^CEBP|^E2F|^HIF|^KLF", .)))



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

ggsave("motif_dotplot_known_sel.pdf", width = 10.5, height = 2.8)
