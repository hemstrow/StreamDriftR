imaf <- runif(10, 0, .5)

#============for the first loci only!==================
N <- 1000

num.minor <- rbinom(1, size = N*2, prob = imaf[1])

nmaf <- num.minor/(2*N)

# loop it for multiple generations!
g <- 100

nmaf <- numeric(g)
nmaf[1] <- imaf[1]

for(i in 2:g){
  nmaf[i] <-  rbinom(1, size = N*2, prob = nmaf[i - 1])/(2*N)
}
plot(nmaf)


#================challenge for next week:=========================
# expand this to do the same process, but with all 10 loci!
# if you have time, try plotting it
# hint: binom(10, 100000, seq(from = .1, to = .5, length.out = 10))/100000



#===============part 2: the part 2ining==========================
# now, track for multiple populations!
# first, use a nested loop to do each pop each generation
# second, vectorize everything to loop only through generations.

# we will need a 3D storage solution(nmaf[gen, locus, population])
# will use an array
npop <- 10
g <- 10

nmaf <- array(NA, c(g, length(imaf), npop))
## access like anything else
nmaf[2,,3] # gen 2, every locus, pop 3
