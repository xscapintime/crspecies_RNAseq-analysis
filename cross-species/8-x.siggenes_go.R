# compare signigicant genes in GO pathway
# between arrested and diapause embyro
# ---------------------------------------

## code from supplementary of
## Scognamiglio, R. et al. (2016)
## ‘Myc Depletion Induces a Pluripotent Dormant State Mimicking Diapause’,
## Cell, 164(4), pp. 668–680. doi: 10.1016/j.cell.2015.12.033.
## -----------------------------------------------------------------------


rm(list = ls())
library(tidyverse)

library(DESeq2)



human_res <- read.table("human_deg.tsv")




library(ggplot2)




data("dxdDiapaused")
dxdDiapaused <- estimateSizeFactors( dxdDiapaused )
filter <- rowSums(counts(dxdDiapaused, normalized=TRUE)) > 0
dxdDiapaused <- dxdDiapaused[filter,]
both <- intersect(rownames(dxdDiapaused), rownames( se ) )

dxdDiapaused <- dxdDiapaused[rownames(dxdDiapaused) %in% both,]
dxdDiapaused <- estimateSizeFactors( dxdDiapaused )
dxdDiapaused <- estimateDispersions(dxdDiapaused)
dxdDiapaused <- nbinomWaldTest( dxdDiapaused )
diapausedRes <- results(dxdDiapaused, independentFiltering = FALSE)

seSub2 <- se[,colData(se)$mix %in% c( "CREplus72h", "ESB8NT0h")]
colData(seSub2) <- droplevels(colData(seSub2))
seSub2 <- seSub2[rownames( seSub2 ) %in% both,]
seSub2 <- estimateSizeFactors( seSub2 )
seSub2 <- estimateDispersions( seSub2 )
seSub2 <- nbinomWaldTest( seSub2 )
resultsdKOvsWT <- results( seSub2, independentFiltering=FALSE )




library(ggplot2)
library(gtable)
library(plyr)
library(dplyr)
#data(mtcars)

scatterBoxplot <- function(dataFrame, highlight=NULL){
if( mean( dataFrame$x[highlight], na.rm=TRUE ) > 0 ){
    col2 <- "#e41a1c80"
}else{
col2 <- "#2166ac80"
}
if( mean( dataFrame$y[highlight], na.rm=TRUE ) > 0 ){
col3 <- "#e41a1c80"
}else{
col3 <- "#2166ac80"
}
dataFrame$col <- ifelse( highlight, col3, "#87878730" )
dataFrame <- dataFrame[rowSums( is.na(dataFrame) ) == 0,]
p1 <- ggplot(dataFrame) +
geom_point(aes(x, y), color="#87878730") +
geom_point(aes(x, y), colour=col3, subset=.(col == col3) ) +
geom_hline(yintercept=0, colour="#00000090", lwd=1.2) +
geom_vline(xintercept=0, colour="#00000090", lwd=1.2) +
xlab(expression(log[2]~"fold change ( Diapaused / E4.5 )")) +
ylab(expression(log[2]~"fold change ( dKO / WT )"))+
theme(plot.margin = unit(c(0.2, 0.2, 0.5, 0.5), "lines"),
panel.background = element_rect(colour="black", fill="white"),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
axis.text=element_text(size=14,color="black"),
axis.title=element_text(size=14,color="black"),
legend.position="none") + ylim(-2.5, 2.5) +
xlim(-7.5, 7.5)
p2 <- ggplot(dataFrame, aes(x = factor(1), y = x)) +
geom_violin(subset=.(col == col3), fill=col2) +
geom_boxplot(outlier.colour = NA, width=.5,subset=.(col == col3),fill=col2) +
geom_hline(yintercept=0, colour="#00000090", lwd=1.2) + ylim(-7.5, 7.5) +
coord_flip() +
theme(axis.text = element_blank(),
axis.title = element_blank(),
axis.ticks = element_blank(),
plot.margin = unit(c(1, 0.2, -0.5, 0.5), "lines"),
panel.background = element_rect(colour="black", fill="white"),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank())
# Vertical marginal boxplot - to appear at the right of the chart
p3 <- ggplot(dataFrame, aes(x = factor(1), y = y)) +
geom_violin(subset=.(col == col3), fill=col3) +
geom_boxplot(outlier.colour = NA, width=.5,subset=.(col == col3), fill=col3) +
geom_hline(yintercept=0, colour="#00000090", lwd=1.2) +
theme(axis.text = element_blank(),
axis.title = element_blank(),
axis.ticks = element_blank(),
plot.margin = unit(c(0.2, 1, 0.5, -0.5), "lines"),
panel.background = element_rect(colour="black", fill="white"),
panel.grid.major = element_blank(),
panel.grid.minor = element_blank()) + ylim(-2.5, 2.5)
gt1 <- ggplot_gtable(ggplot_build(p1))
gt2 <- ggplot_gtable(ggplot_build(p2))
gt3 <- ggplot_gtable(ggplot_build(p3))
maxWidth <- unit.pmax(gt1$widths[2:3], gt2$widths[2:3])
maxHeight <- unit.pmax(gt1$heights[4:5], gt3$heights[4:5])
# Set the maximums in the gtables for gt1, gt2 and gt3
gt1$widths[2:3] <- as.list(maxWidth)
gt2$widths[2:3] <- as.list(maxWidth)
gt1$heights[4:5] <- as.list(maxHeight)
gt3$heights[4:5] <- as.list(maxHeight)
gt <- gtable(widths = unit(c(7, 1.5), "null"),
height = unit(c(1.5, 7), "null"))
gt <- gtable_add_grob(gt, gt1, 2, 1)
gt <- gtable_add_grob(gt, gt2, 1, 1)
gt <- gtable_add_grob(gt, gt3, 2, 2)
# And render the plot
grid.newpage()
}
