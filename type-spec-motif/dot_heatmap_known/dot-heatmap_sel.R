# dot heatmap
# looks like @malikPluripotencyReprogrammingCompetent2019 Fig 2b
# ---------------------------------------------------------------

rm(list = ls())

library(tidyverse)

### load homer tables
load("known_motif_all.Rdata")

dat <- dat_kn

#### ==================================================== ####
#### fig in the article
# FOX, ATF, MYC, MYB, RUNX, HNF, PPAR, CEBP, E2F, HIF and KLF

dat_sel <- dat %>% filter(across(fam, ~ grepl("FOX|ATF|MYC|MYB|RUNX|HNF|PPAR|CEBP|E2F|HIF|^KLF", .)))


# remove double motif
dat_sel <- dat_sel %>% filter(across(name, ~ !grepl(":", .)))


## plot ==================================================================================
## colors
mycols <- c("#4B80DB", "#dc1a91","#dc1a84", "#dc1a64", "#DC1E1A")
morecols <- (grDevices::colorRampPalette(mycols))(125) ### expand the platte!


# selected
theme_set(
    cowplot::theme_minimal_grid(line_size = 0.3) +
    theme(axis.line  = element_line(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1 , size = 7, margin = margin(t = 0, b = 1)),
        axis.text.y = element_text(size = 7),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.key.size = unit(.2, "cm"),
        plot.margin = margin(t = 1.5, r = .5, b = .5, l = .5, unit = "cm")
    )
)

# pick by hand...
sel <- c("KLF3", "KLF14", "E2F4", "E2F", "ATF1", "ATF4", "CEBP", "MYB", "c-MYC", "FOXH1", "FOXO1", "HIF-1B", "HIF-1A", "HNF4A", "RUNX", "PPARA")


dat_sel %>% mutate(`P-value (-log10)` = -log10(P.value)) %>% filter(ortholog %in% sel) %>%
    ggplot(aes(x = ortholog, y = type, fill = `P-value (-log10)`, size = percent)) +
    geom_point(alpha = 1, shape = 21, colour = "black", stroke = .3) +
    ylab("") + xlab("") +
    scale_fill_gradientn(colours = morecols) +
    #scale_size(trans = "reverse") +
    scale_y_discrete(limits = rev,
        labels = rev(c("Type I DE", "Type II DE", "Type III DE"))) +
    labs(size = "% Motif in Target", fill = bquote(-log["10"] (P-value))) +
    cowplot::panel_border(color  = "black", size = .6) + coord_fixed()
    # scale_size(range = c(0, 5))


ggsave("motif_dotplot_known_sel.pdf", width = 5.1, height = 2.1)
