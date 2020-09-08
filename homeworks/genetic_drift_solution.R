imaf <- runif(1000, 0, .5)

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


N <- 100
g <- 1000
nmaf <- matrix(NA, g, length(imaf))
nmaf[1,] <- imaf
for(i in 2:g){
  for(j in 1:length(imaf)){
    nmaf[i,j] <- rbinom(1, N*2, prob = nmaf[i-1,j])/(N*2)
  }
}

N <- 100
nmaf <- matrix(NA, g, length(imaf))
nmaf[1,] <- imaf
for(i in 2:g){
  nmaf[i,] <- rbinom(length(imaf), N*2, prob = nmaf[i-1,])/(N*2)
}
