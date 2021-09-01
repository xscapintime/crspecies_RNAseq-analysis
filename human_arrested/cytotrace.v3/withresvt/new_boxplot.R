# boxplot for subsetted cytotrace values
# --------------------------------------

rm(list = ls())

library(tidyverse)


## cytotrace table
cyt <- read.table("CytoTRACE_plot_table.txt", header = T, row.names = 1, sep = "\t")

## subset
sub <- cyt %>% filter(!Phenotype %in% c("E5", "E6","E7", "Blastocyst", "Late blastocyst"))
sub$Phenotype <- factor(sub$Phenotype, levels = c("Arrested Type I", "Zygote", "2C", "4C", "Oocyte",
                        "Arrested Type II", "8C", "E3", "Arrested Type III", "Morula", "Arrested Type R", "E4"))


## draw boxplot
library(ggplot2)
#library(hrbrthemes)


cols <- c("#9E0142", "#C0274A", "#DC494C", "#F06744", "#F88D51", "#FDB466", "#FDD380", "#FEDC56",
        "#FDD700", "#EFE95E", "#D7EF9B", "#B2E0A2", "#88CFA4", "#5FBAA8", "#3F96B7", "#4272B2", "#5E4FA2")


theme_set(
    theme_test() +
    theme(
        legend.position = "none",
        plot.title = element_text(size = 11)
    ))


p <- ggplot(sub, aes(x = Phenotype, y = CytoTRACE, color = Phenotype, fill = Phenotype))
p +  geom_boxplot(alpha = 0.4) +
    geom_jitter(size = 0.3, width = 0.25) +
    ylab("Predicted ordering by CytoTRACE") + xlab("") +
    scale_color_manual(values = adjustcolor(cols, alpha.f = 0.25)[1:12], aesthetics = "fill") +
    scale_color_manual(values = cols[seq_along(unique(sub$Phenotype))], aesthetics = "color") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggsave(width = 7.6, height = 7, filename = "cytotrace_boxplot_withres_part.pdf")