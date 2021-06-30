system.time(
        tmp1 <- t(sapply(split(human_exp_mat[, -29], human_exp_mat$hgnc_symbol),
        function(x) sapply(x, mean)))
)

system.time(
        tmp2 <- aggregate( . ~ hgnc_symbol, data = human_exp_mat[, -29], mean)
)


system.time(
        tmp3 <- human_exp_mat[, -29]
        %>% group_by(hgnc_symbol)
        %>% summarise_all(mean)
)

any(tmp1[,1:28] != tmp2[, 2:29])
any(tmp1[,1:28] != tmp3[, 2:29])
any(tmp2 != tmp3)