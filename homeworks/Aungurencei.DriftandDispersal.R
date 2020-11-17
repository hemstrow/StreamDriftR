
dm <- mc2d::rmultinomial(nrow(out), N, t(out))
dm 
N2 <- colSums(dm)
dm <- mc2d::rmultinomial(nrow(out), N2, t(out))


# hints:
tm <- matrix(1:9, 3, 3) # dispersal numbers, rows are source pop
afs <- matrix(runif(9, 0, .5), 3, 3) # afs, rows are loci
draws <- rbinom(3*3*3, rep(tm*2, each = 3), afs) # this will flip coins correctly
# for locus 1, destination 1, the new allele frequency is:
sum(draws[c(1, 4, 7)])/ colSums(tm)[1]
