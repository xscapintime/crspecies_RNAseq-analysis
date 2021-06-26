library(visdat)


png("figs/human_navis.png")
vis_dat(human_res)
dev.off()


png("figs/mouse_navis.png")
vis_dat(mouse_res)
dev.off()