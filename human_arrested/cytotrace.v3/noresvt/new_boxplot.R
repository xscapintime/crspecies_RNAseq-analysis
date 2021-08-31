# boxplot for subsetted cytotrace values
# --------------------------------------

rm(list = ls())

library(tidyverse)


## cytotrace table
cyt <- read.table("CytoTRACE_plot_table.txt", header = T, row.names = 1, sep = "\t")

## subset
sub <- cyt %>% filter(!Phenotype %in% c("E4", "E5", "E6","E7", "Blastocyst", "Late blastocyst"))
sub$Phenotype <- factor(sub$Phenotype, levels = c("Oocyte", "Zygote", "2C", "Arrested Type I", "4C",
                        "Arrested Type II", "8C", "E3", "Arrested Type III", "Morula"))


## draw boxplot
library(ggplot2)
library(hrbrthemes)


cols <- c("#9E0142", "#C2294A", "#DF4D4B", "#F46D43", "#FA9856", "#FDBE6E", "#FEE08B", "#FEDA2E",
            "#F6E132", "#E6F497", "#BEE5A0", "#94D4A4", "#66C2A5", "#439BB5", "#4075B4", "#5E4FA2")


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
    scale_color_manual(values = adjustcolor(cols, alpha.f = 0.25)[1:10], aesthetics = "fill") +
    scale_color_manual(values = cols[1:10], aesthetics = "color") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggsave(width = 7, height = 6, filename = "cytotrace_boxplot_nores_part.pdf")





## cytotrace way
#---------------
pheno_levels <-as.character(unique(sub$Phenotype))

pdf("CytoTRACE_Boxplot_new.pdf", width = length(pheno_levels) * 1.5, height = 6)

par(mar = c(log(max(nchar(pheno_levels)) + 2, 2) * 3.25, 6, 2, 1), xpd = NA)

boxplot(sub$CytoTRACE ~ sub$Phenotype, outline = F,
        las = 2, xaxt = "n", yaxt = "n", staplelwd = 2,
        medlwd = 2, whisklty = 1, border = cols, col = adjustcolor(cols,
        alpha.f = 0.25), ylab = "Predicted ordering by CytoTRACE",
        xlab = "", cex.lab = 1.75, frame.plot = F,
        xlim = c(0,length(pheno_levels) + 0.5), ylim = c(min(sub$CytoTRACE), 1))

mtext("Cell phenotypes", cex = 1.75, side = 1, line = round(max(nchar(pheno_levels)/2, 1) + 1))

axis(1, pos = 0, at = 1:length(pheno_levels), labels = FALSE)
text(x = seq_along(pheno_levels), y = -0.05, srt = 45, adj = 1, labels = pheno_levels, xpd = TRUE, cex = 1.5)

segments(0, 0, length(pheno_levels) + 0.5, 0)
axis(2, pos = 0, at = seq(0, 1, 0.2), labels = format(seq(0, 1, 0.2), nsmall = 1), las = 2, cex.axis = 1.5)
segments(0, 0, 0, 1)

stripchart(sub$CytoTRACE ~ sub$Phenotype, vertical = TRUE,
            method = "jitter", add = TRUE, pch = 16, col = cols,
            bg = "white", lwd = 1, cex = 2/round(log(length(sub$CytoTRACE), 10)))

dev.off()
