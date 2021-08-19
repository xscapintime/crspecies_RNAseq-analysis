# position score as corrdinate
# p-val as dot size?
# -----------------------------


rm(list = ls())

library(tidyverse)


# position score
load("position_score.Rdata")


# BH p-val
load("manova_bh_adjpval.Rdata")



## plot data
dat <- cbind(human_ps, mouse_ps, bh_adjpval) %>% data.frame

dat$group <- ifelse(dat$mouse_ps * dat$human_ps > 0,
                    ifelse(abs(dat$mouse_ps) >= 0.1 & abs(dat$human_ps) >= 0.1, "corr", "other"),
                    ifelse(abs(dat$mouse_ps) >= 0.1 & abs(dat$human_ps) >= 0.1, "anti", "other"))


dat$sig <- ifelse(dat$bh_adjpval < 0.05, "sig", "no")


# get pathway category
# bpidx <- row.names(dat) %>% grep(pattern = "GOBP_")
# mfidx <- row.names(dat) %>% grep(pattern = "GOMF_")
# ccidx <- row.names(dat) %>% grep(pattern = "GOCC_")
# kgidx <- row.names(dat) %>% grep(pattern = "KEGG_")
# diaidx <- row.names(dat) %>% grep(pattern = "Diapause_")

dat$cat <- unlist(lapply(X = row.names(dat), FUN = function(x) {
                return(strsplit(x, split = "_")[[1]][1])}))

write.table(dat, file = "manova-position.txt", quote = F, sep = "\t")


# correlation test
sigs <- dat %>% filter(bh_adjpval < 0.03) %>% select(human_ps, mouse_ps)

cor.test(sigs$mouse_ps, sigs$human_ps, method = "pearson")
cor.test(sigs$mouse_ps, sigs$human_ps, method = "spearman")



## dotplot
library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(ggthemes)
library(ggrepel)


##########
# for all
#########

theme_set(theme_few() +
    theme(legend.justification = c(1, 0),
        legend.position = c(0.98, 0.01),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.key.size = unit(18, "pt"))
)



#formula <- y ~ x

p <- ggplot(dat %>% filter(bh_adjpval < 0.02),
            aes(x = mouse_ps, y = human_ps))
p + geom_point(alpha = .6, aes(color = group, size = bh_adjpval,
                shape = cat)) + #shape = sig

    scale_size(trans = "reverse") +
    scale_color_manual(values = c("#D43F3A", "#5CB85C", "#7C878E")) +

    # geom_text(data = dat %>%
    #           filter(group == "corr" &
    #           padj.x <= 0.01 & padj.y <= 0.01),
    #           mapping = aes(label = pathway), size = 3) +

    labs(x = "Diapause / E4.5", y = "Arrested / 8C") +
    xlim(-1, 1) + ylim(-1, 1) +
    geom_vline(xintercept = 0, linetype = "dashed", size = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", size = 0.6) +
    guides(color = F, shape = guide_legend(title = NULL),
            size = guide_legend(title = "q-value")) +
    scale_shape_discrete(breaks = c("GOBP", "GOMF", "GOCC", "KEGG", "Diapause")) +
    stat_cor(data = dat %>% filter(bh_adjpval < 0.02) %>% select(human_ps, mouse_ps),
                method = "spearman")
    #geom_smooth(method = "lm")

ggsave(width = 7.6, height = 7.6, filename = "../figs/2d_allpath.png")




### Add text

p <- ggplot(dat %>% filter(bh_adjpval < 0.02),
            aes(x = mouse_ps, y = human_ps))
p + geom_point(alpha = .6, aes(color = group, size = bh_adjpval,
                shape = cat)) + #shape = sig

    scale_size(trans = "reverse") +
    scale_color_manual(values = c("#D43F3A", "#5CB85C", "#7C878E")) +

    geom_text_repel(data = dat %>%
                filter(abs(dat$mouse_ps * dat$human_ps) > 0.125 & group == "corr"),
                mapping = aes(label = dat %>%
                                filter(abs(dat$mouse_ps * dat$human_ps) > 0.125 & group == "corr") %>%
                                row.names), size = 2.5, segment.size = 0.2, segment.color = "#DCDCDC") +

    labs(x = "Diapause / E4.5", y = "Arrested / 8C") +
    xlim(-1, 1) + ylim(-1, 1) +
    geom_vline(xintercept = 0, linetype = "dashed", size = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", size = 0.6) +
    guides(color = F, shape = guide_legend(title = NULL),
            size = guide_legend(title = "q-value")) +
    scale_shape_discrete(breaks = c("GOBP", "GOMF", "GOCC", "KEGG", "Diapause")) +
    stat_cor(data = dat %>% filter(bh_adjpval < 0.02) %>% select(human_ps, mouse_ps),
                method = "spearman")
    #geom_smooth(method = "lm")

ggsave(width = 7.6, height = 7.6, filename = "../figs/2d_allpath_corr.png")



p <- ggplot(dat %>% filter(bh_adjpval < 0.02),
            aes(x = mouse_ps, y = human_ps))
p + geom_point(alpha = .6, aes(color = group, size = bh_adjpval,
                shape = cat)) + #shape = sig

    scale_size(trans = "reverse") +
    scale_color_manual(values = c("#D43F3A", "#5CB85C", "#7C878E")) +

    geom_text_repel(data = dat %>%
                filter(abs(dat$mouse_ps * dat$human_ps) > 0.125 & group == "anti"),
                mapping = aes(label = dat %>%
                                filter(abs(dat$mouse_ps * dat$human_ps) > 0.125 & group == "anti") %>%
                                row.names), size = 2.5, segment.size = 0.4, segment.color = "#DCDCDC") +

    labs(x = "Diapause / E4.5", y = "Arrested / 8C") +
    xlim(-1, 1) + ylim(-1, 1) +
    geom_vline(xintercept = 0, linetype = "dashed", size = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", size = 0.6) +
    guides(color = F, shape = guide_legend(title = NULL),
            size = guide_legend(title = "q-value")) +
    scale_shape_discrete(breaks = c("GOBP", "GOMF", "GOCC", "KEGG", "Diapause")) +
    stat_cor(data = dat %>% filter(bh_adjpval < 0.02) %>% select(human_ps, mouse_ps),
                method = "spearman")
    #geom_smooth(method = "lm")

ggsave(width = 7.6, height = 7.6, filename = "../figs/2d_allpath_anti.pdf")




p <- ggplot(dat %>% filter(bh_adjpval < 0.02),
            aes(x = mouse_ps, y = human_ps))
p + geom_point(alpha = .6, aes(color = group, size = bh_adjpval,
                shape = cat)) + #shape = sig

    scale_size(trans = "reverse") +
    scale_color_manual(values = c("#D43F3A", "#5CB85C", "#7C878E")) +

    geom_text_repel(data = dat %>%
                filter(abs(dat$mouse_ps * dat$human_ps) > 0.18),
                mapping = aes(label = dat %>%
                                filter(abs(dat$mouse_ps * dat$human_ps) > 0.18) %>%
                                row.names), size = 2.5, segment.size = 0.5, segment.color = "#DCDCDC") +

    labs(x = "Diapause / E4.5", y = "Arrested / 8C") +
    xlim(-1, 1) + ylim(-1, 1) +
    geom_vline(xintercept = 0, linetype = "dashed", size = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", size = 0.6) +
    guides(color = F, shape = guide_legend(title = NULL),
            size = guide_legend(title = "q-value")) +
    scale_shape_discrete(breaks = c("GOBP", "GOMF", "GOCC", "KEGG", "Diapause")) +
    stat_cor(data = dat %>% filter(bh_adjpval < 0.02) %>% select(human_ps, mouse_ps),
                method = "spearman")
    #geom_smooth(method = "lm")

ggsave(width = 7.6, height = 7.6, filename = "../figs/2d_allpath_seltext.pdf")




##########################
# for 4 diapause gene sets
##########################

dat2 <- dat %>% filter(cat == "Diapause")
dat2$exp <- c("up", "dn", "up", "dn")



p <- ggplot(dat2, aes(x = mouse_ps, y = human_ps))
p + geom_point(alpha = .6, aes(size = bh_adjpval, color = exp)) +

    scale_size(trans = "reverse") +
    scale_color_manual(values = c("#186eb4", "#ec4848")) +

    geom_text_repel(aes(label = row.names(dat2))) +

    labs(x = "Diapause / E4.5", y = "Arrested / 8C") +
    xlim(-1, 1) + ylim(-1, 1) +
    geom_vline(xintercept = 0, linetype = "dashed", size = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", size = 0.6) +
    guides(color = F, shape = guide_legend(title = NULL),
            size = guide_legend(title = "q-value"))
    # stat_cor(data = dat2 %>% select(human_ps, mouse_ps),
    #             method = "spearman")
    #geom_smooth(method = "lm")

ggsave(width = 7.6, height = 7.6, filename = "../figs/2d_diapause.png")

