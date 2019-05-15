library(ggplot2)


#==================growth and dispersal==================
L<-100 #100 meters
n<- 30 #number of spots
m<- 1 #mean location they disperse to
a<- 2 #width of the dispersal 


#spread and growth in one generation

r<-1 #growth rate
K<-3000 #carrying capacity
##BH model 10 generations
gen <- 20



# run the function
output <- the.function(L, n, m, a, r, K, gen)


# plot
output <- as.data.frame(output)
output$gen <- 1:gen

m.o <- reshape2::melt(output, id.vars = "gen") # converts to long for. The id.vars are the variables that identify the quality of the rows that we are intersted in.
colnames(m.o)[2:3] <- c("position", "N") # give new column names
m.o$gen <- as.factor(m.o$gen) # change the classes of the columns so that ggplot doesn't throw a fit. Generally, things that are truely categories should be factors. Things that are truely continuous should be numeric.
m.o$position <- as.numeric(as.character(m.o$position)) # to convert factor to numeric, need to go through a character intermediate. It's screwy.
ggplot(m.o, aes(x = position, y = N, color = gen)) + geom_line() + theme_bw()
# ggplot general format:
# data, then your aesthetics (x, y, color, fill, shape, ect)--the things that define your plot
# then, add on another function with + that says what kind of plot you want.



#==================genetic drift======================
# next goal: do genetic drift at one location.
# example data:
as <- readRDS("drift_example_data.RDS") # snp dataset

# assume these are all at the same location. Simulate genetic drift for gen generations.
# assume no growth, constant pop size of N.

# example for one row, one generation! The allele identity (A/C/G/T) is not important, don't track it. Just track p and q across generations.
ex.row <- as[3,]
N <- 420

# calculate p, the frequency of the most common allele
p <- matrixStats::rowMaxs(as)/rowSums(as)

# here is how you draw alleles for one site!
p1 <- p[3]
np.next <- rbinom(N, 2, p1) # how many p alleles each individal in the next gen got.
np.n.total <- sum(np.next)
# possibly usefull alternative
np.n.total <- rbinom(1, 2*N, p1) # since we don't care about individal genotype, this is like mixing them all together.
npf <- np.n.total/(2*N)

# now do this for every row all at once. Avoid a loop here. (rbinom is vectorized for the probability!) rbinom(2,1, c(.001, .9))
# Then, do this for multiple generations. Need a loop here.
# make a single function that does this. Inputs should be only the parameters.