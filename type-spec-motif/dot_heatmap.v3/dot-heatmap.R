# dot heatmap
# looks like @malikPluripotencyReprogrammingCompetent2019 Fig 2b
# ---------------------------------------------------------------

rm(list = ls())

library(tidyverse)

## load homer tables

files <- list.files("../filtering.v2", pattern = "homer-\\w*.txt", recursive = F)
homer <- lapply(files, function(file) { read.csv(paste0("../filtering.v2/", file), sep = "\t", header = T)})
names(homer) <- unlist(lapply(X = files, FUN = function(x) {return(strsplit(x, split = ".txt")[[1]][1])}))


## known
typespec_kn <- lapply(homer[-(1:2)], function(x) select(x, `Motif.Name`, `P.value.known`, `percent.known`) %>% na.omit)

typespec_kn <- Map(function(dat, grp) cbind(dat, group = grp), typespec_kn,
                unlist(lapply(X = names(typespec_kn), FUN = function(x) {return(strsplit(x, split = "homer-")[[1]][2])})))

dat_kn <- Reduce(rbind, typespec_kn)
dat_kn$group <- factor(dat_kn$group, levels = c("typeI_up", "typeII_up", "typeIII_up", "typeI_dn", "typeII_dn", "typeIII_dn"))
dat_kn$tf <- unlist(lapply(X = dat_kn$Motif.Name, FUN = function(x) {return(strsplit(x, split = "/")[[1]][1])}))


## de novo
typespec_novo <- lapply(homer[-(1:2)], function(x) filter(x, is.na(P.value.known)) %>% select(`Motif.Name`, `P.value.novo`, `percent.novo`, 'score'))
typespec_novo <- typespec_novo[sapply(typespec_novo, nrow) > 0]

typespec_novo <- Map(function(dat, grp) cbind(dat, group = grp), typespec_novo,
                unlist(lapply(X = names(typespec_novo), FUN = function(x) {return(strsplit(x, split = "homer-")[[1]][2])})))

dat_novo <- Reduce(rbind, typespec_novo)
dat_novo$group <- factor(dat_novo$group, levels = c("typeI_up", "typeII_up", "typeIII_up", "typeI_dn", "typeII_dn", "typeIII_dn"))
dat_novo$tf <- unlist(lapply(X = dat_novo$Motif.Name, FUN = function(x) {return(strsplit(x, split = "/")[[1]][1])}))


# combine
colnames(dat_kn) <- c("Motif.Name", "P.value", "percent", "group", "tf")
colnames(dat_novo) <- c("Motif.Name", "P.value", "percent", "score", "group", "tf")

dat <- rbind(dat_kn, dat_novo[-4])


## names
dat$name <- unlist(lapply(X = dat$tf, FUN = function(x) {return(strsplit(x, split = "\\(")[[1]][1])}))
dat$name[grepl("_", dat$name)] <- unlist(lapply(X = dat$name[grepl("_", dat$name)], FUN = function(x) {return(strsplit(x, split = "_")[[1]][2])}))
#dat <- dat %>% group_by(group, name) %>% summarise(mean(P.value), mean(percent)) # kinda wrong


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

# oder by pval and percentage
tf_order <- (dat_sel %>% group_by(fam, group) %>% arrange(desc(percent), P.value, by_group = T))$name %>% unique()
dat_sel$name <- factor(dat_sel$name, levels = tf_order)





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
    ggplot(aes(x = name, y = group, fill = `P-value (-log10)`, size = percent)) +
    geom_point(alpha = 1, shape = 21, colour = "black", stroke = .5) +
    ylab("") + xlab("") +
    scale_fill_gradientn(colours = morecols) +
    #scale_size(trans = "reverse") +
    scale_y_discrete(limits = rev,
        labels = rev(c("Type I up", "Type II up", "Type III up", "Type I down", "Type II down", "Type III down"))) +
    labs(size = "% Motif in Target", fill = bquote(-log["10"] (P-value))) +
    cowplot::panel_border(color  = "black", size = .8)

ggsave("motif_dotplot_sel.pdf", width = 10.5, height = 3)





# all type, all motif

dat$de <- unlist(lapply(X = as.character(dat$group), FUN = function(x) {return(strsplit(x, split = "_")[[1]][2])}))

# oder by pval and percentage
tf_order <- (dat %>% group_by(de, group) %>% arrange(P.value, desc(percent), by_group = T))$name %>% unique()
dat$name <- factor(dat$name, levels = tf_order)



theme_set(
    cowplot::theme_minimal_grid() +
    theme(axis.line  = element_line(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1 , size = 5),
        axis.text.y = element_text(size = 6),
        #axis.ticks = element_blank(),
        panel.border = element_blank(),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        legend.key.size = unit(.3, "cm"),
        plot.margin = margin(t = 1.5, r = 1, b = .5, l = 1, unit = "cm")
    )
)

dat <- dat %>% group_by(group, name) %>% filter(percent == max(percent))

dat %>% mutate(`P-value (-log10)` = -log10(P.value)) %>% #filter(percent > 10)

    ggplot(aes(x = name, y = group, fill = `P-value (-log10)`, size = percent)) +
    geom_point(alpha = 1, shape = 21, colour = "black", stroke = .3) +
    ylab("") + xlab("") +
    scale_fill_gradientn(colours = morecols) +
    scale_y_discrete(limits = rev,
        labels = rev(c("Type I up", "Type II up", "Type III up", "Type I down", "Type II down", "Type III down"))) +
    labs(size = "% Motif in Target", fill = bquote(-log["10"] (P-value))) +
    cowplot::panel_border(color  = "black", size = .8) +
    scale_size(range = c(0, 3))


ggsave("motif_dotplot_alltype.pdf", width = 20, height = 2.5)
