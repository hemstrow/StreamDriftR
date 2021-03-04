#================challenge for next week:=========================
# expand this to do the same process, but with all 10 loci!
# if you have time, try plotting it
# hint: rbinom(10, 100000, seq(from = .1, to = .5, length.out = 10))/100000

nmaf <- matrix(NA, g, length(imaf))
nmaf[1, ] <- imaf
for (i in 2:g){
  nmaf[i,] <- rbinom(length(imaf), size = N*2, prob = nmaf[i-1,])/(N *2)
}

matplot(nmaf)