
dm <- mc2d::rmultinomial(nrow(out), N, t(out))
dm 
N2 <- colSums(dm)
dm <- mc2d::rmultinomial(nrow(out), N2, t(out))

