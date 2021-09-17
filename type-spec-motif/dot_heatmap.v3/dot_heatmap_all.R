# dot heatmap
# looks like @malikPluripotencyReprogrammingCompetent2019 Fig 2b
# ---------------------------------------------------------------

rm(list = ls())

library(tidyverse)

### load homer tables

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
# typespec_novo <- lapply(homer[-(1:2)], function(x) filter(x, is.na(P.value.known)) %>% select(`Motif.Name`, `P.value.novo`, `percent.novo`, 'score'))
# typespec_novo <- typespec_novo[sapply(typespec_novo, nrow) > 0]

# typespec_novo <- Map(function(dat, grp) cbind(dat, group = grp), typespec_novo,
#                 unlist(lapply(X = names(typespec_novo), FUN = function(x) {return(strsplit(x, split = "homer-")[[1]][2])})))

# dat_novo <- Reduce(rbind, typespec_novo)
# dat_novo$group <- factor(dat_novo$group, levels = c("typeI_up", "typeII_up", "typeIII_up", "typeI_dn", "typeII_dn", "typeIII_dn"))
# dat_novo$tf <- unlist(lapply(X = dat_novo$Motif.Name, FUN = function(x) {return(strsplit(x, split = "/")[[1]][1])}))


## combine
colnames(dat_kn) <- c("Motif.Name", "P.value", "percent", "group", "tf")
# colnames(dat_novo) <- c("Motif.Name", "P.value", "percent", "score", "group", "tf")

## remove the score column
# dat <- rbind(dat_kn, dat_novo[-4])
dat <- dat_kn


## names
dat$name <- unlist(lapply(X = dat$tf, FUN = function(x) {return(strsplit(x, split = "\\(")[[1]][1])}))
# dat$name[grepl("_", dat$name)] <- unlist(lapply(X = dat$name[grepl("_", dat$name)], FUN = function(x) {return(strsplit(x, split = "_")[[1]][2])}))
#dat <- dat %>% group_by(group, name) %>% summarise(mean(P.value), mean(percent)) # kinda wrong

## group and order

# if same TF and same group, take the one with higher percentage
dat <- dat %>% group_by(group, name) %>% filter(P.value == min(P.value))

# Tf family
reg <- c("MYB","^AP-1|^Atf|^ATF|^BATF","^Bach","^bHLHE","^CEBP","^CTCF","^E2F","^Egr","^ELF|^Elf","^Elk","^ETS|^Ets|^ETV|^Etv","^EWS","^FOX|^Fox","^Fra","^Gata|^GATA","^GFY","^HIF","^Hox|^HOX","^IRF\\d","^Jun","^KLF|^Klf","^Mef","^NFY","^NFkB","^Nkx","^NPAS","^NRF","^PHOX2A|^Phox2a","^PRDM","^Rfx|^RFX","^Sox","^Sp\\d|^SP\\d","^Spi|PU.1$","^Stat|^STAT","^TEAD","^USF|^Usf","^ZBTB","^ZNF|^Znf|^ZSCAN")
fam <- c("MYB","AP-1/ATF","BACH","bHLHE","CEBP","CTCF","E2F","EGR","ELF","ELK","ETS","EWS","FOX","FRA","GATA","GFY","HIF","HOX","IRF","JUN","KLF","MEF","NFY","NFKB","NKX","NPAS","NRF","PHOX2A","PRDM","RFX","SOX","SP","SPI","STAT","TEAD","USF","ZBTB","ZNF")

dat$fam <- dat$name
for (i in 1:length(reg)) {
    dat$fam[grepl(reg[i], dat$name)] <- fam[i]
}

# split 1 and 2,3
dat$type <- unlist(lapply(X = as.character(dat$group), FUN = function(x) {return(strsplit(x, split = "_")[[1]][1])}))
dat$typegroup <- ifelse(dat$type == "typeI", "type1", "type23")

# de
dat$de <- unlist(lapply(X = as.character(dat$group), FUN = function(x) {return(strsplit(x, split = "_")[[1]][2])}))

# make TF family factors
fam_order <- (dat %>% group_by(typegroup, fam) %>% summarise(mean(percent), min(P.value))
    %>% arrange(`min(P.value)`, desc(`mean(percent)`), .by_group = T))$fam %>% unique()
dat$fam <- factor(dat$fam, levels = fam_order)

# change the name to human orthologs
# not compelete
dat$ortholog <- toupper(dat$name)

# oder by pval and percentage
tf_order <- (dat %>% group_by(typegroup, fam) %>% arrange(P.value, desc(percent), .by_group = T))$ortholog %>% unique()
dat$ortholog <- factor(dat$ortholog, levels = tf_order)


## save the table
save(dat, file = "motif_known_all.Rdata")
write.table(dat, file = "motif_knwon_all.txt", quote = F, sep = "\t", row.names = F, col.names = T)


# oder by pval and percentage
# tf_order <- (dat %>% group_by(de, group) %>% arrange(P.value, desc(percent), by_group = T))$name %>% unique()
# dat$name <- factor(dat$name, levels = tf_order)



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
