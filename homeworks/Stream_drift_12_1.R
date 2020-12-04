imaf <- matrix(runif(200, 0, .5), 10, 20) 

npop <- 20
g <- 15
ipops <- seq(from= 30, to= 1000, length.out = 20)
ipops <- floor(ipops)

source("R/stream_drift.R")
out <- A(1, -5, 50, 99)
N <- numeric(99)
N[1] <- 100

ipops <- mc2d::rmultinomial(nrow(out), N, t(out)) #dispersal matrix

nmaf <- array(NA, c(g, nrow(imaf), npop))#setting parameter 
nmaf[1,,] <-imaf

for (i in 2:g){
  ipops <- mc2d::rmultinomial(nrow(out), N, t(out)) #dispersal matrix
  rbinom(n= nrow(imaf)* npop, size= ipops, prob= imaf) # how many alleles that fish is dispersing at each locus
  DAT <- data.frame(
    dest = rep(c(1:npop), each=nrow(imaf)),
    locus = rep(c(1:nrow(imaf)), npop),
    allele_count = rbinom(n= nrow(imaf)* npop, size= ipops, prob= imaf),
    N = rep(ipops, each= nrow(imaf)))
  New_allele_freq <- 
    tapply(DAT$allele_count,list(DAT$dest,DAT$locus),  FUN=sum)/
    tapply(DAT$N *2, list(DAT$dest,DAT$locus),  FUN=sum)
  nmaf[i,,] <- rbinom(n= nrow(imaf)* npop, size= ipops, prob= New_allele_freq)
}

