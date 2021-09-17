# dot heatmap
# looks like @malikPluripotencyReprogrammingCompetent2019 Fig 2b
# ---------------------------------------------------------------

rm(list = ls())

library(tidyverse)

### load homer tables
load("known_motif_all.Rdata")


### ==========================================================================================
## plot
## colors
mycols <- c("#4B80DB", "#dc1a91","#dc1a84", "#dc1a64", "#DC1E1A")
morecols <- (grDevices::colorRampPalette(mycols))(125) ### expand the platte!


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


dat %>% mutate(`P-value (-log10)` = -log10(P.value)) %>% #filter(percent > 10)

    ggplot(aes(x = ortholog, y = group, fill = `P-value (-log10)`, size = percent)) +
    geom_point(alpha = 1, shape = 21, colour = "black", stroke = .3) +
    ylab("") + xlab("") +
    scale_fill_gradientn(colours = morecols) +
    scale_y_discrete(limits = rev,
        labels = rev(c("Type I up", "Type II up", "Type III up", "Type I down", "Type II down", "Type III down"))) +
    labs(size = "% Motif in Target", fill = bquote(-log["10"] (P-value))) +
    cowplot::panel_border(color  = "black", size = .8) +
    scale_size(range = c(0, 2.5))


ggsave("motif_dotplot_knwon_all.pdf", width = 18, height = 2.8)
