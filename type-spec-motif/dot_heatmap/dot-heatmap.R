# dot heatmap
# looks like @malikPluripotencyReprogrammingCompetent2019 Fig 2b
# ---------------------------------------------------------------

rm(list = ls())

library(tidyverse)

## load homer tables

files <- list.files("../filtering", pattern = "homer-\\w*.txt", recursive = F)
homer <- lapply(files, function(file) { read.csv(paste0("../filtering/", file), sep = "\t", header = T)})
names(homer) <- unlist(lapply(X = files, FUN = function(x) {return(strsplit(x, split = ".txt")[[1]][1])}))


## maybe no de novo result
typespec <- lapply(homer[-(1:2)], function(x) select(x, `Motif.Name`, `P.value.known`, `percent.known`) %>% na.omit)

typespec <- lapply(typespec, function(x) x %>% mutate(group = names(x)))

typespec <- Map(function(dat, grp) cbind(dat, group = grp), typespec,
                unlist(lapply(X = names(typespec), FUN = function(x) {return(strsplit(x, split = "homer-")[[1]][2])})))

dat <- Reduce(rbind, typespec)
dat$group <- factor(dat$group, levels = c("typeI_up", "typeII_up", "typeIII_up", "typeI_dn", "typeII_dn", "typeIII_dn"))
dat$tf <- unlist(lapply(X = dat$Motif.Name, FUN = function(x) {return(strsplit(x, split = "/")[[1]][1])}))


## colors
mycols <- c("#4B80DB","#dc1a84", "#dc1a64", "#DC1E1A")
morecols <- (grDevices::colorRampPalette(mycols))(50) ### expand the platte!


# all
dat %>% mutate(`P-value (-log10)` = -log10(P.value.known)) %>%
    ggplot(aes(x = tf, y = group, fill = `P-value (-log10)`, size = percent.known)) +
    geom_point(alpha = 1, shape = 21, colour = "black") +
    cowplot::theme_minimal_grid() +
    theme(axis.line  = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1 , size = 7)) +
    ylab("") + xlab("") +
    theme(axis.ticks = element_blank()) +
    theme(, legend.title = element_text(size = 10), legend.text = element_text(size = 9)) +
    scale_fill_gradientn(colours = morecols) +
    #scale_size(trans = "reverse") +
    scale_y_discrete(limits = rev,
        labels = rev(c("Type I up", "Type II up", "Type III up", "Type I down", "Type II down", "Type III down"))) +
    guides(size = guide_legend(title = "% Motif in Target"))

ggsave("motif_dotplot_all.pdf", width = 15, height = 7)


# type I only
dat %>% mutate(`P-value (-log10)` = -log10(P.value.known)) %>% filter(group %in% c("typeI_up", "typeI_dn")) %>%
    ggplot(aes(x = tf, y = group, fill = `P-value (-log10)`, size = percent.known)) +
    geom_point(alpha = 1, shape = 21, colour = "black") +
    cowplot::theme_minimal_grid() +
    theme(axis.line  = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1 , size = 7)) +
    ylab("") + xlab("") +
    theme(axis.ticks = element_blank()) +
    theme(, legend.title = element_text(size = 10), legend.text = element_text(size = 9)) +
    scale_fill_gradientn(colours = morecols) +
    #scale_size(trans = "reverse") +
    scale_y_discrete(limits = rev,
        labels = rev(c("Type I up", "Type I down"))) +
    guides(size = guide_legend(title = "% Motif in Target"))

ggsave("motif_dotplot_typeI.pdf", width = 9, height = 5)


# type II only
dat %>% mutate(`P-value (-log10)` = -log10(P.value.known)) %>% filter(group %in% c("typeII_up", "typeII_dn")) %>%
    ggplot(aes(x = tf, y = group, fill = `P-value (-log10)`, size = percent.known)) +
    geom_point(alpha = 1, shape = 21, colour = "black") +
    cowplot::theme_minimal_grid() +
    theme(axis.line  = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1 , size = 7)) +
    ylab("") + xlab("") +
    theme(axis.ticks = element_blank()) +
    theme(, legend.title = element_text(size = 10), legend.text = element_text(size = 9)) +
    scale_fill_gradientn(colours = morecols) +
    #scale_size(trans = "reverse") +
    scale_y_discrete(limits = rev,
        labels = rev(c("Type II up", "Type II down"))) +
    guides(size = guide_legend(title = "% Motif in Target"))

ggsave("motif_dotplot_typeII.pdf", width = 9, height = 5)


# type III only
dat %>% mutate(`P-value (-log10)` = -log10(P.value.known)) %>% filter(group %in% c("typeIII_up", "typeIII_dn")) %>%
    ggplot(aes(x = tf, y = group, fill = `P-value (-log10)`, size = percent.known)) +
    geom_point(alpha = 1, shape = 21, colour = "black") +
    cowplot::theme_minimal_grid() +
    theme(axis.line  = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1 , size = 7)) +
    ylab("") + xlab("") +
    theme(axis.ticks = element_blank()) +
    theme(, legend.title = element_text(size = 10), legend.text = element_text(size = 9)) +
    scale_fill_gradientn(colours = morecols) +
    #scale_size(trans = "reverse") +
    scale_y_discrete(limits = rev,
        labels = rev(c("Type III up", "Type III down"))) +
    guides(size = guide_legend(title = "% Motif in Target"))

ggsave("motif_dotplot_typeIII.pdf", width = 9, height = 5)



# up only
dat %>% mutate(`P-value (-log10)` = -log10(P.value.known)) %>% filter(across(group, ~ !grepl("up", .))) %>%
    ggplot(aes(x = tf, y = group, fill = `P-value (-log10)`, size = percent.known)) +
    geom_point(alpha = 1, shape = 21, colour = "black") +
    cowplot::theme_minimal_grid() +
    theme(axis.line  = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1 , size = 7)) +
    ylab("") + xlab("") +
    theme(axis.ticks = element_blank()) +
    theme(, legend.title = element_text(size = 10), legend.text = element_text(size = 9)) +
    scale_fill_gradientn(colours = morecols) +
    #scale_size(trans = "reverse") +
    scale_y_discrete(limits = rev,
        labels = rev(c("Type I up", "Type II up", "Type III up"))) +
    guides(size = guide_legend(title = "% Motif in Target"))

ggsave("motif_dotplot_up.pdf", width = 9, height = 5)


# down only
dat %>% mutate(`P-value (-log10)` = -log10(P.value.known)) %>% filter(across(group, ~ !grepl("dn", .))) %>%
    ggplot(aes(x = tf, y = group, fill = `P-value (-log10)`, size = percent.known)) +
    geom_point(alpha = 1, shape = 21, colour = "black") +
    cowplot::theme_minimal_grid() +
    theme(axis.line  = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1 , size = 7)) +
    ylab("") + xlab("") +
    theme(axis.ticks = element_blank()) +
    theme(, legend.title = element_text(size = 10), legend.text = element_text(size = 9)) +
    scale_fill_gradientn(colours = morecols) +
    #scale_size(trans = "reverse") +
    scale_y_discrete(limits = rev,
        labels = rev(c("Type I down", "Type II down", "Type III down"))) +
    guides(size = guide_legend(title = "% Motif in Target"))

ggsave("motif_dotplot_dn.pdf", width = 9, height = 5)
