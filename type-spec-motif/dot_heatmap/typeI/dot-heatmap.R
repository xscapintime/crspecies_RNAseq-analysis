gene_cluster %>% filter(Gene %in% markers) %>%
    mutate(`% Expressing` = (cell_exp_ct/cell_ct) * 100) %>%
    filter(count > 0, `% Expressing` > 1) %>%
    ggplot(aes(x=cluster, y = Gene, color = count, size = `% Expressing`)) +
    geom_point() +
    cowplot::theme_cowplot() +
    theme(axis.line  = element_blank()) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ylab('') +
    theme(axis.ticks = element_blank()) +
    scale_color_gradientn(colours = viridis::viridis(20), limits = c(0,4), oob = scales::squish, name = 'log2 (count + 1)')