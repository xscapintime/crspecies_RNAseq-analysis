# 2D plot based on path FC
# position socre from 9-1
# GO BP for a look
# ------------------------


rm(list = rm())

library(tidyverse)


human_bp <- read.csv("human_bp_fc_ps.txt", sep = "\t", header = T)
mouse_bp <- read.csv("mouse_bp_fc_ps.txt", sep = "\t", header = T)


dat <- inner_join(mouse_bp, human_bp, by = "GO.ID") %>%
                    dplyr::select(GO.ID, Term.x, ps.x, ps.y,
                                    classicFisher.x, classicFisher.y)



dat$group <- ifelse(dat$ps.x * dat$ps.y > 0,
                    ifelse(abs(dat$ps.x) >= 0.8 & abs(dat$ps.y) >= 0.8, "corr", "other"),
                    ifelse(abs(dat$ps.x) >= 0.8 & abs(dat$ps.y) >= 0.8, "anti", "other"))


dat$sig <- ifelse(dat$classicFisher.x < 0.01 & dat$classicFisher.y < 0.01, "bothsig",
                    ifelse(dat$classicFisher.x < 0.01, "msig",
                            ifelse(dat$classicFisher.y < 0.01, "hsig", "no")))



### dot plot 
## x: moue path position score
## y: human path position score

theme_set(
  theme_classic() +
    theme(legend.position = "right"))


p <- ggplot(dat, #%>% filter(padj.x <= 0.05 | padj.y <= 0.05),
            aes(x = ps.x, y = ps.y))
p + geom_point(alpha = .6, size = 1, aes(color = group, shape = sig)) +
                                #size = dat$size.x/dat$size.y)) + # size are the same

    scale_color_manual(values = c("#D43F3A", "#5CB85C", "#7C878E")) +

    # geom_text(data = dat %>%
    #           filter(group == "corr" &
    #           padj.x <= 0.01 & padj.y <= 0.01),
    #           mapping = aes(label = pathway), size = 3) +

    labs(title = "2D GOBP pathways",
          x = "Arrested / 8C pathway FC", y = "Diapause / E4.5 pathway FC")

ggsave(width = 7.6, height = 7.6, filename = "figs/2dgobp_pathfcps.png")
