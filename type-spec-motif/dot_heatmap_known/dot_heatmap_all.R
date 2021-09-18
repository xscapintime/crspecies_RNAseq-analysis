# dot heatmap
# looks like @malikPluripotencyReprogrammingCompetent2019 Fig 2b
# ---------------------------------------------------------------

rm(list = ls())

library(tidyverse)

### load homer tables
load("known_motif_all.Rdata")

dat <- dat_kn

### ==========================================================================================
## plot
## colors
mycols <- c("#4B80DB", "#dc1a91","#dc1a84", "#dc1a64", "#DC1E1A")
morecols <- (grDevices::colorRampPalette(mycols))(125) ### expand the platte!


theme_set(
    cowplot::theme_minimal_grid(line_size = 0.3) +
    theme(axis.line  = element_line(),
        axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1 , size = 4, margin = margin(t = 0, b = 1)),
        axis.text.y = element_text(size = 5),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 5),
        legend.key.size = unit(.2, "cm"),
        plot.margin = margin(t = 1.5, r = 1, b = 1, l = 1, unit = "cm")
    )
)

# to wrap the plot
# options(download.file.method = "wininet")
# devtools::install_github("wilkox/ggwrap")
# library(ggwrap)


dat %>% mutate(`P-value (-log10)` = -log10(P.value)) %>% #filter(percent > 10)
        ggplot(aes(x = ortholog, y = type, fill = `P-value (-log10)`, size = percent)) +
        geom_point(alpha = 1, shape = 21, colour = "black", stroke = .3) +
        ylab("") + xlab("") +
        scale_fill_gradientn(colours = morecols) +
        scale_y_discrete(limits = rev,
            labels = rev(c("Type I DE", "Type II DE", "Type III DE"))) +
        labs(size = "% Motif in Target", fill = bquote(-log["10"] (P-value))) +
        cowplot::panel_border(color  = "black", size = .5) +
        scale_size(range = c(0, 2.6))

ggsave("motif_dotplot_knwon_all.pdf", width = 20, height = 2.1)


# wp <- ggwrap(p, 3)
# ggsave("motif_dotplot_knwon_all_wrap.pdf", plot = wp, width = 8, height = 3.5)


## divide to 3?4? rows
d <- dat %>% mutate(`P-value (-log10)` = -log10(P.value))

d1 <- d %>% filter(ortholog %in% levels(dat$ortholog)[1:91])
d2 <- d %>% filter(ortholog %in% levels(dat$ortholog)[92:182])
d3 <- d %>% filter(ortholog %in% levels(dat$ortholog)[183:271])



theme_set(
    cowplot::theme_minimal_grid(line_size = 0.3) +
    theme(axis.line  = element_line(),
        axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1 , size = 4, margin = margin(t = 0, b = 1)),
        axis.text.y = element_text(size = 5),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 5),
        legend.key.size = unit(.2, "cm"),
        plot.margin = margin(t = 0, r = 1, b = 1, l = 0, unit = "cm")
    )
)


# library(ggh4x)

p1 <- d1 %>% ggplot(aes(x = ortholog, y = type, fill = `P-value (-log10)`, size = percent)) +
        geom_point(alpha = 1, shape = 21, colour = "black", stroke = .3) +
        ylab("") + xlab("") +
        scale_fill_gradientn(colours = morecols, limits = c(min(d$`P-value (-log10)`), max(d$`P-value (-log10)`))) +
        scale_size_continuous(limits = c(min(d$percent), max(d$percent))) +
        scale_y_discrete(limits = rev,
            labels = rev(c("Type I DE", "Type II DE", "Type III DE"))) +
        labs(size = "% Motif in Target", fill = bquote(-log["10"] (P-value))) +
        cowplot::panel_border(color  = "black", size = .6) +
        # force_panelsizes(rows = 1, cols = 10) +
        coord_fixed() +
        scale_size(range = c(0, 3))

p2 <- d2 %>% ggplot(aes(x = ortholog, y = type, fill = `P-value (-log10)`, size = percent)) +
        geom_point(alpha = 1, shape = 21, colour = "black", stroke = .3) +
        ylab("") + xlab("") +
        scale_fill_gradientn(colours = morecols, limits = c(min(d$`P-value (-log10)`), max(d$`P-value (-log10)`))) +
        scale_size_continuous(limits = c(min(d$percent), max(d$percent))) +
        scale_y_discrete(limits = rev,
            labels = rev(c("Type I DE", "Type II DE", "Type III DE"))) +
        labs(size = "% Motif in Target", fill = bquote(-log["10"] (P-value))) +
        cowplot::panel_border(color  = "black", size = .6) +
        # force_panelsizes(rows = 1, cols = 10) +
        coord_fixed() +
        scale_size(range = c(0, 3))

p3 <- d3 %>% ggplot(aes(x = ortholog, y = type, fill = `P-value (-log10)`, size = percent)) +
        geom_point(alpha = 1, shape = 21, colour = "black", stroke = .3) +
        ylab("") + xlab("") +
        scale_fill_gradientn(colours = morecols, limits = c(min(d$`P-value (-log10)`), max(d$`P-value (-log10)`))) +
        scale_size_continuous(limits = c(min(d$percent), max(d$percent))) +
        scale_y_discrete(limits = rev,
            labels = rev(c("Type I DE", "Type II DE", "Type III DE"))) +
        labs(size = "% Motif in Target", fill = bquote(-log["10"] (P-value))) +
        cowplot::panel_border(color  = "black", size = .6) +
        # force_panelsizes(rows = 1, cols = 10) +
        coord_fixed() +
        scale_size(range = c(0, 3))


gridp <- cowplot::plot_grid(p1, p2, p3, nrow = 3, align = "v")
ggsave("motif_grid_all.pdf", gridp)
