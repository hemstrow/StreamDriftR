imaf <- runif(20*10, 0, .5)
imaf <- matrix(imaf, 10, 20) # rows are loci, columsn are populations

g <- 15 # fifteen generations

nmaf <- array(NA, c(g, nrow(imaf), ncol(imaf)))
## access like anything else
nmaf[2,,3] # gen 2, every locus, pop 3
nmaf[1,,] <- imaf

nmaf[1,,]

# here's the goal:
# for one generation, use the dispersal code you wrote a while back to figure out how many individuals move from each location
# to each other location. Then, figure the new allele frequencies at each location after these individuals all move around. Each
# individual carries to allele at each loci to the new populations, but draws them from their source populations.
# use a single rbinom call, make sure to align the number of individuals (argument 2) in each group (d1s1, d1s2, etc) to the allele frequencies 
# (argument 2, d1s1l1, d1s1l2, ..., d20s20l10)
