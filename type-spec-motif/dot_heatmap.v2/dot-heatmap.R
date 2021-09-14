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


# combne
colnames(dat_kn) <- c("Motif.Name", "P.value", "percent", "group", "tf")
colnames(dat_novo) <- c("Motif.Name", "P.value", "percent", "score", "group", "tf")

dat <- rbind(dat_kn, dat_novo[-4])

## colors
mycols <- c("#4B80DB","#dc1a84", "#dc1a64", "#DC1E1A")
morecols <- (grDevices::colorRampPalette(mycols))(50) ### expand the platte!


# all
dat %>% mutate(`P-value (-log10)` = -log10(P.value)) %>% filter(percent > 10 & `P-value (-log10)` > 3) %>%
    ggplot(aes(x = tf, y = group, fill = `P-value (-log10)`, size = percent)) +
    geom_point(alpha = 1, shape = 21, colour = "black") +
    cowplot::theme_minimal_grid() +
    theme(axis.line  = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1 , size = 6)) +
    ylab("") + xlab("") +
    theme(axis.ticks = element_blank()) +
    theme(, legend.title = element_text(size = 10), legend.text = element_text(size = 9)) +
    scale_fill_gradientn(colours = morecols) +
    #scale_size(trans = "reverse") +
    scale_y_discrete(limits = rev,
        labels = rev(c("Type I up", "Type II up", "Type III up", "Type I down", "Type II down", "Type III down"))) +
    guides(size = guide_legend(title = "% Motif in Target"))

ggsave("motif_dotplot_all.pdf", width = 20, height = 7)


## =============================
## need collapse into clusters
dat <- dat %>% mutate(`P-value (-log10)` = -log10(P.value)) %>% filter(percent > 10 & `P-value (-log10)` >= 3)

myb <- dat %>% filter(across(Motif.Name, ~ grepl("Myb|MYB", .))) %>% group_by(group) %>% summarise(mean(P.value), mean(percent)) %>% mutate(clusname = "MYB")
atf <- dat %>% filter(across(Motif.Name, ~ grepl("Atf|ATF", .))) %>% group_by(group) %>% summarise(mean(P.value), mean(percent)) %>% mutate(clusname = "ATF")
bach <- dat %>% filter(across(Motif.Name, ~ grepl("Bach|BACH", .))) %>% group_by(group) %>% summarise(mean(P.value), mean(percent)) %>% mutate(clusname = "BACH")
bhlhe <- dat %>% filter(across(Motif.Name, ~ grepl("bHLHE", .))) %>% group_by(group) %>% summarise(mean(P.value), mean(percent)) %>% mutate(clusname = "bHLH")
cebp <- dat %>% filter(across(Motif.Name, ~ grepl("CEBP", .))) %>% group_by(group) %>% summarise(mean(P.value), mean(percent)) %>% mutate(clusname = "CEBP")
e2f <- dat %>% filter(across(Motif.Name, ~ grepl("E2F", .))) %>% group_by(group) %>% summarise(mean(P.value), mean(percent)) %>% mutate(clusname = "E2F")
egr <- dat %>% filter(across(Motif.Name, ~ grepl("Egr", .))) %>% group_by(group) %>% summarise(mean(P.value), mean(percent)) %>% mutate(clusname = "EGR")
#elf <- dat %>% filter(across(Motif.Name, ~ grepl("Elf|ELF", .))) %>% group_by(group) %>% summarise(mean(P.value), mean(percent))
#elk <- dat %>% filter(across(Motif.Name, ~ grepl("Elk|ELK", .))) %>% group_by(group) %>% summarise(mean(P.value), mean(percent))
ets <- dat %>% filter(across(Motif.Name, ~ grepl("ETS|ETV", .))) %>% group_by(group) %>% summarise(mean(P.value), mean(percent)) %>% mutate(clusname = "ETS")
fox <- dat %>% filter(across(Motif.Name, ~ grepl("Forkhead", .))) %>% group_by(group) %>% summarise(mean(P.value), mean(percent)) %>% mutate(clusname = "FOX")
fra <- dat %>% filter(across(Motif.Name, ~ grepl("Fra", .))) %>% group_by(group) %>% summarise(mean(P.value), mean(percent)) %>% mutate(clusname = "FRA")
gata <- dat %>% filter(across(Motif.Name, ~ grepl("Gata|GATA", .))) %>% group_by(group) %>% summarise(mean(P.value), mean(percent)) %>% mutate(clusname = "GATA")
hif <- dat %>% filter(across(Motif.Name, ~ grepl("HIF", .))) %>% group_by(group) %>% summarise(mean(P.value), mean(percent)) %>% mutate(clusname = "HIF")
hox <- dat %>% filter(across(Motif.Name, ~ grepl("^Hox|^HOX", .))) %>% group_by(group) %>% summarise(mean(P.value), mean(percent)) %>% mutate(clusname = "HOX")
irf <- dat %>% filter(across(Motif.Name, ~ grepl("^IRF\\d", .))) %>% group_by(group) %>% summarise(mean(P.value), mean(percent)) %>% mutate(clusname = "IRF")
klf <- dat %>% filter(across(Motif.Name, ~ grepl("^KLF|^Klf", .))) %>% group_by(group) %>% summarise(mean(P.value), mean(percent)) %>% mutate(clusname = "KLF")
nkx <- dat %>% filter(across(Motif.Name, ~ grepl("Nkx", .))) %>% group_by(group) %>% summarise(mean(P.value), mean(percent)) %>% mutate(clusname = "NKX")
nrf <- dat %>% filter(across(Motif.Name, ~ grepl("NRF", .))) %>% group_by(group) %>% summarise(mean(P.value), mean(percent)) %>% mutate(clusname = "NRF")
phox2a <- dat %>% filter(across(Motif.Name, ~ grepl("PHOX2A|Phox2a", .))) %>% group_by(group) %>% summarise(mean(P.value), mean(percent)) %>% mutate(clusname = "PHOX2A")
prdm <- dat %>% filter(across(Motif.Name, ~ grepl("PRDM", .))) %>% group_by(group) %>% summarise(mean(P.value), mean(percent)) %>% mutate(clusname = "PRDM")
sox <- dat %>% filter(across(Motif.Name, ~ grepl("Sox", .))) %>% group_by(group) %>% summarise(mean(P.value), mean(percent)) %>% mutate(clusname = "SOX")
sp <- dat %>% filter(across(Motif.Name, ~ grepl("Sp\\d", .))) %>% group_by(group) %>% summarise(mean(P.value), mean(percent)) %>% mutate(clusname = "SP")
stat <- dat %>% filter(across(Motif.Name, ~ grepl("Stat|STAT", .))) %>% group_by(group) %>% summarise(mean(P.value), mean(percent)) %>% mutate(clusname = "STAT")
tbox <- dat %>% filter(across(Motif.Name, ~ grepl("T-box", .))) %>% group_by(group) %>% summarise(mean(P.value), mean(percent)) %>% mutate(clusname = "T-box")
tead <- dat %>% filter(across(Motif.Name, ~ grepl("TEAD", .))) %>% group_by(group) %>% summarise(mean(P.value), mean(percent)) %>% mutate(clusname = "TEAD")
usf <- dat %>% filter(across(Motif.Name, ~ grepl("USF|Usf", .))) %>% group_by(group) %>% summarise(mean(P.value), mean(percent)) %>% mutate(clusname = "USF")
zbtb <- dat %>% filter(across(Motif.Name, ~ grepl("ZBTB", .))) %>% group_by(group) %>% summarise(mean(P.value), mean(percent)) %>% mutate(clusname = "ZBTB")
znf <- dat %>% filter(across(Motif.Name, ~ grepl("ZNF", .))) %>% group_by(group) %>% summarise(mean(P.value), mean(percent)) %>% mutate(clusname = "ZNF")


# rest of it
clps_idx <- c("Myb|MYB|Atf|ATF|Bach|BACH|bHLHE|CEBP|E2F|Egr|ETS|ETV|Forkhead|Fra|Gata|GATA|HIF|^Hox|^HOX|^IRF\\d|^KLF|^Klf|Nkx|NRF|PHOX2A|Phox2a|PRDM|Sox|Sp\\d|Stat|STAT|T-box|TEAD|USF|Usf|ZBTB|ZNF")

dat_rest <- dat %>% filter(across(Motif.Name, ~!grepl(clps_idx, .)))

dat_rest$name <- unlist(lapply(X = dat_rest$tf, FUN = function(x) {return(strsplit(x, split = "\\(")[[1]][1])}))
dat_rest$name[grepl("_", dat_rest$name)] <- unlist(lapply(X = dat_rest$name[grepl("_", dat_rest$name)], FUN = function(x) {return(strsplit(x, split = "_")[[1]][2])}))
dat_rest <- dat_rest %>% group_by(group, name) %>% summarise(mean(P.value), mean(percent))

dat_clps <- rbind(myb,atf,bach,bhlhe,cebp,e2f,egr,ets,fox,fra,gata,hif,hox,irf,klf,nkx,nrf,phox2a,prdm,sox,sp,stat,tbox,tead,usf,zbtb,znf)

colnames(dat_rest) <- c("group", "name", "P.value", "percent")
colnames(dat_clps) <- c("group", "P.value", "percent", "name")
new_dat <- rbind(dat_clps, dat_rest)



## ===============================
# all type
new_dat %>% mutate(`P-value (-log10)` = -log10(P.value)) %>% filter(percent > 10 & `P-value (-log10)` >= 3) %>%
    ggplot(aes(x = name, y = group, fill = `P-value (-log10)`, size = percent)) +
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

ggsave("motif_dotplot_alltype.pdf", width = 12, height = 7)




# type I only
new_dat %>% mutate(`P-value (-log10)` = -log10(P.value)) %>% filter(group %in% c("typeI_up", "typeI_dn")) %>%
    ggplot(aes(x = name, y = group, fill = `P-value (-log10)`, size = percent)) +
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
new_dat %>% mutate(`P-value (-log10)` = -log10(P.value)) %>% filter(group %in% c("typeII_up", "typeII_dn")) %>%
    ggplot(aes(x = name, y = group, fill = `P-value (-log10)`, size = percent)) +
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
new_dat %>% mutate(`P-value (-log10)` = -log10(P.value)) %>% filter(group %in% c("typeIII_up", "typeIII_dn")) %>%
    ggplot(aes(x = name, y = group, fill = `P-value (-log10)`, size = percent)) +
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
new_dat %>% mutate(`P-value (-log10)` = -log10(P.value)) %>% filter(across(group, ~ !grepl("up", .))) %>%
    ggplot(aes(x = name, y = group, fill = `P-value (-log10)`, size = percent)) +
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

ggsave("motif_dotplot_up.pdf", width = 10, height = 5)


# down only
new_dat %>% mutate(`P-value (-log10)` = -log10(P.value)) %>% filter(across(group, ~ !grepl("dn", .))) %>%
    ggplot(aes(x = name, y = group, fill = `P-value (-log10)`, size = percent)) +
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

ggsave("motif_dotplot_dn.pdf", width = 10, height = 5)



# =============================
## common and specific
new_dat$de <- unlist(lapply(X = as.character(new_dat$group), FUN = function(x) {return(strsplit(x, split = "_")[[1]][2])}))
new_dat$type <- unlist(lapply(X = as.character(new_dat$group), FUN = function(x) {return(strsplit(x, split = "_")[[1]][1])}))


# common up for 3 types
up_common_idx <- Reduce(intersect,list((new_dat %>% filter(de == "up" & type == "typeI"))$name,
                    (new_dat %>% filter(de == "up" & type == "typeII"))$name,
                    (new_dat %>% filter(de == "up" & type == "typeIII"))$name))

# common dn for 3 types
dn_common_idx <- Reduce(intersect,list((new_dat %>% filter(de == "dn" & type == "typeI"))$name,
                    (new_dat %>% filter(de == "dn" & type == "typeII"))$name,
                    (new_dat %>% filter(de == "dn" & type == "typeIII"))$name))


# common up

new_dat %>% mutate(`P-value (-log10)` = -log10(P.value)) %>% filter(name %in% up_common_idx & de == "up") %>%
    ggplot(aes(x = name, y = group, fill = `P-value (-log10)`, size = percent)) +
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

ggsave("motif_dotplot_comm_up.pdf", width = 7, height = 5)



# common dn

new_dat %>% mutate(`P-value (-log10)` = -log10(P.value)) %>% filter(name %in% dn_common_idx & de == "dn") %>%
    ggplot(aes(x = name, y = group, fill = `P-value (-log10)`, size = percent)) +
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

ggsave("motif_dotplot_comm_dn.pdf", width = 7, height = 5)



# type 1 specifc

new_dat %>% mutate(`P-value (-log10)` = -log10(P.value)) %>%
    filter(!name %in% unique(c(up_common_idx, dn_common_idx)) & type == "typeI") %>%
    ggplot(aes(x = name, y = group, fill = `P-value (-log10)`, size = percent)) +
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

ggsave("motif_dotplot_typeI_specifc.pdf", width = 7, height = 5)


# type 2 specifc

new_dat %>% mutate(`P-value (-log10)` = -log10(P.value)) %>%
    filter(!name %in% unique(c(up_common_idx, dn_common_idx)) & type == "typeII") %>%
    ggplot(aes(x = name, y = group, fill = `P-value (-log10)`, size = percent)) +
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

ggsave("motif_dotplot_typeII_specifc.pdf", width = 7, height = 5)


# type 3 specifc

new_dat %>% mutate(`P-value (-log10)` = -log10(P.value)) %>%
    filter(!name %in% unique(c(up_common_idx, dn_common_idx)) & type == "typeIII") %>%
    ggplot(aes(x = name, y = group, fill = `P-value (-log10)`, size = percent)) +
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

ggsave("motif_dotplot_typeIII_specifc.pdf", width = 7, height = 5)



## similar?
typeIII_up_id <- (new_dat %>% mutate(`P-value (-log10)` = -log10(P.value)) %>%
    filter(!name %in% unique(c(up_common_idx, dn_common_idx)) & group == "typeIII_up"))$name

typeII_up_id <- (new_dat %>% mutate(`P-value (-log10)` = -log10(P.value)) %>%
    filter(!name %in% unique(c(up_common_idx, dn_common_idx)) & group == "typeII_up"))$name

typeI_up_id <- (new_dat %>% mutate(`P-value (-log10)` = -log10(P.value)) %>%
    filter(!name %in% unique(c(up_common_idx, dn_common_idx)) & group == "typeI_up"))$name


# of course zero
Reduce(intersect, list(typeI_up_id, typeII_up_id, typeIII_up_id))

# typw 1 vs type 2, 32 vs 22
Reduce(intersect, list(typeI_up_id, typeII_up_id))
# [1] "HIF"       "NRF"       "STAT"      "USF"       "ZBTB"
# [6] "Arnt:Ahr"  "CLOCK"     "COUP-TFII" "EKLF"      "NFYA"
# [11] "Rbpj1"     "WT1"

# type 2 vs type 3, 22 vs 26
Reduce(intersect, list(typeII_up_id, typeIII_up_id))
# [1] "SOX"       "c-Jun-CRE" "Mef2b"     "Pit1+1bp"

# type 1 vs type 3, 32 vs 26
Reduce(intersect, list(typeI_up_id, typeIII_up_id))
# [1] "EGR"   "EAR2"  "GSC"   "HNF4a" "JunB"  "Otx2"



# down
typeIII_dn_id <- (new_dat %>% mutate(`P-value (-log10)` = -log10(P.value)) %>%
    filter(!name %in% unique(c(up_common_idx, dn_common_idx)) & group == "typeIII_dn"))$name

typeII_dn_id <- (new_dat %>% mutate(`P-value (-log10)` = -log10(P.value)) %>%
    filter(!name %in% unique(c(up_common_idx, dn_common_idx)) & group == "typeII_dn"))$name

typeI_dn_id <- (new_dat %>% mutate(`P-value (-log10)` = -log10(P.value)) %>%
    filter(!name %in% unique(c(up_common_idx, dn_common_idx)) & group == "typeI_dn"))$name


# of course zero
Reduce(intersect, list(typeI_dn_id, typeII_dn_id, typeIII_dn_id))

# typw 1 vs type 2, 19 vs 13
Reduce(intersect, list(typeI_dn_id, typeII_dn_id))
# [1] "CEBP"  "HIF"   "NRF"   "STAT"  "HLF"   "MITF"  "Rbpj1"

# type 2 vs type 3, 13 vs 9
Reduce(intersect, list(typeII_dn_id, typeIII_dn_id))
# [1] "Prop1"

# type 1 vs type 3, 19 vs 9
Reduce(intersect, list(typeI_dn_id, typeIII_dn_id))
# [1] "EGR"   "CLOCK"