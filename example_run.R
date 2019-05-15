library(ggplot2); library(reshape2); source("R/stream_drift.R"); source("R/sn2.R")
# parameters
r <- 2 #set r
k <- 3000 #set K
Tf <- 50 #set number of gens
l <- 10 #set number of loci


# get an edge df. Can be gotten from GIStoEdge!
bs <- data.frame(branch.id = c("001", "002", "003", "004", "005", "006"), 
                 length = c(100, 25, 100, 75, 25, 20),
                 start.node = c("C", "B", "A", "D", "E", "D"),
                 start.weight = c(1, .5, 1, 1, 1, 1),
                 end.node = c("A", "A", "F", "B", "B", "C"),
                 end.weight = c(1,1,.5,.2,.8,1), stringsAsFactors = FALSE)


# convert to a dispersal matrix
dm <- Edge.to.Matrix(bs, 2, 5, 1)

# normalize so that all individuals end up somewhere (it's a close system). Can change latter
for (i in 1:ncol(dm[[1]])){
  dm[[1]][,i] <- dm[[1]][,i]/sum(dm[[1]][,i])
}

# set the initial pop sizes. As given, this assumes that individuals start as upstream as possible in branch
n0 <- numeric(nrow(dm[[2]]))
n0[301] <- 100


# grab the positions for the matrix entries
xs <- dm[[2]]$xs

# run the stream drift function! currently producing NaNs, need to work on it.
out <- stream.drift(n0, l, NULL, r, NULL, NULL, k, Tf, in.mat = dm[[1]], in.xs = xs)

# add branch metadata back in
out$dat$branch <- dm[[2]]$branch

# melt and plot
mdat <- melt(out$dat, id.vars = c("N", "xs", "branch"))
colnames(mdat)[c(4,5)] <- c("loci", "maf")
ggplot(mdat, aes(x = xs, y = maf, color = loci)) + geom_line() + theme_bw() + facet_wrap(~branch)
