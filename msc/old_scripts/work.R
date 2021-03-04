l <- 10
maf <- .2
mafsd <- .2
n <- 100

lmaf <- rnorm(n = l, mean = maf, sd = mafsd) #generate mafs from mafsd and mean maf
lmaf <- lmaf[lmaf > 0 & lmaf < 1] #truncate
while(length(lmaf) < l){ #while lmaf stays too short...
  lmaf <- c(lmaf, rnorm(n = l, mean = maf, sd = mafsd)) #add more draws
  lmaf <- lmaf[lmaf > 0 & lmaf < 1] #truncate
}
lmaf <- lmaf[1:l]

gen_genos <- function(tmaf){
  q <- tmaf
  draws <- rbinom(n, 2, q)
  draws_geno <- ifelse(draws == 0,"0202", ifelse(draws == 1, "0102", "0101"))
  return(draws_geno)
}

test <- lapply(FUN = gen_genos, lmaf)

((length(grep(pattern = "0101", test[[1]]))*2) + length(grep(pattern = "0102", test[[1]])))/200
v <- list()
for(i in 1:length(test)){
  v <- c(v, ((length(grep(pattern = "0101", test[[i]]))*2) + length(grep(pattern = "0102", test[[i]])))/200)
}
mean(unlist(v))

test2 <- as.data.frame(matrix(unlist(test), nrow = l, ncol = n, byrow = TRUE))

lapply(X = "ind", FUN = paste0, 1:30)



test <- make_pop(n = 10, l = 100)
v <- list()
for(i in 1:nrow(test)){
  v <- c(v, ((length(grep(pattern = "0101", test[i,]))*2) + length(grep(pattern = "0102", test[i,])))/20)
}


n <- 1000

lmaf <- make_pop(n,l)
lmaf <- rep(0.5, length = l)

#make_pop(n, l, test, pop.out = TRUE)


tf <- 1000
N <- matrix(NA, tf, l)
N[1,] <- lmaf

for(i in 1:(tf - 1)){
  N[i+1,] <- make_pop(n,l,N[i,])
}

matplot(N, type = "l")







#Working on the influence matrix and such
a <- 2 #set a
R <- 4 #set R
m <- 1 #set m
L <- 10 #set L
n <- 100 #set n
l <- 10 #set number of loci
lmaf <- rep(0.5, length = l) #get random mafs from uniform dist

#get the input matrix
L10 <- A(a, R, m, L, n, "both")
L10xs <- unlist(L10[[2]])
L10mat <- unlist(L10[[1]])


#iterate through time modeling pop growth and dispersal with the matrix
Tf <- 10
N <- matrix(NA, ncol = length(L10xs), nrow = Tf)
N[1,] <- ifelse(L10xs < -9 & L10xs > -10, 1, 0)

for(i in 1:(Tf-1)){
  N[i+1,] <- L10mat%*%N[i,]
}


##plot
library(ggplot2)
#melt, move iNo df for easy plotting
library(reshape2)
N <- t(N)
N <- as.data.frame(N)
colnames(N) <- c(1:Tf) #set column names
N[,"x"] <- L10xs #set a column of positions
N_m <- melt(N, id.vars = c("x")) #melt data
colnames(N_m) <- c("x", "t", "value")

plot_1 <- ggplot(data = N_m, aes(x = x, y = value, color = t)) + geom_line() +
  xlab("Spatial Position") + ylab("N")


#iterate through maf, need to iterate through pop sizes as this happens
Maf1 <- matrix(NA, nrow = length(L10xs), ncol= Tf)
Maf1[,1] <- lmaf[1]

N <- ifelse(L10xs < -9 & L10xs > -10, 10, 0) 

for(i in 1:(Tf-1)){
  N <- L10mat%*%N
  #cat("i is:", i, "\n")
  for(j in 1:nrow(Maf1)){
    #("j is:", j, "N[j]:", N[j], "maf:", Maf1[j,i], "i+1:", i+1, "\n")
    Maf1[j,i+1] <- make_pop(N[j], 1, maf = Maf1[j,i], pop.out = FALSE)
    #print(Maf1[j,i+1])
  }
}

test_maf <- c(rep(.45, 10), rep(.1, 90))



#How do I weight mafs?
#cells are the influence of column on row
#sum influence on position one is row 1
#need to weight maf impact by all spots on row
#need relative impact first
rel_weights <- L10mat[1,]/sum(L10mat[1,]) #get relative weights
#get weighted maf from each
w_mafs <- test_maf*rel_weights

#sum weighted maf
weight1 <- sum(w_mafs)

#how do I do this for more rows/locations at once?
#start with a function to weight
w.maf <- function(weights, mafs){
  out <- sum(mafs*weights/sum(weights))
  return(out)
}
w.maf(L10mat[1,], test_maf)
#weights change, mafs don't
test <- apply(L10mat, 1, w.maf, mafs = test_maf) 
#works, this takes each row of the influence matrix A and the mafs for each position corresponding
#to the rows of that matrix and calculates the weighted incoming maf for those locations.

#do I get the same results if I use an influence matrix with no growth? I should.
#When I removed the R part of the A function, I got exactly the same thing! No problems here.

#so, the full function for calculating mafs is, with w.maf not redefined...
make_wmafs <- function(in.mafs, i.mat){
  out <- apply(i.mat, 1, w.maf, mafs = in.mafs)
}

test2 <- make_wmafs(test_maf, L10mat)


##Ok, let's try drifting the a constant population according to these weights through a few generations.
# To get a starting point, allow a few gens of dispersal, then hold those constant and drift.
Tf <- 10
N <- matrix(NA, ncol = length(L10xs), nrow = Tf)
N[1,] <- ifelse(L10xs < -9 & L10xs > -10, 1, 0)

for(i in 1:(Tf-1)){
  N[i+1,] <- L10mat%*%N[i,]
}

N <- N[7,]

Tf <- 100
#with no drift, just spread
L_m <- matrix(NA, nrow = length(N), ncol = Tf)
L_m[,1] <- test_maf
for(i in 1:(Tf-1)){
  this_maf <- L_m[,i]
  this_maf <- make_wmafs(in.mafs = this_maf, i.mat = L10mat)
  L_m[,i+1] <- this_maf
}
matplot(t(L_m), type = "l")

#with drift
L_m <- matrix(NA, nrow = length(N), ncol = Tf)
L_m[,1] <- test_maf
for(i in 1:(Tf-1)){
  this_maf <- L_m[,i]
  this_maf <- make_wmafs(in.mafs = this_maf, i.mat = L10mat)
  #cat("Weighted_maf:")
  #print(this_maf)
  for(j in 1:length(this_maf)){#drift in each pop
    this_maf[j] <- make_pop(N[j], 1, this_maf[j])
  }
  #cat("After drift:")
  #print(this_maf)
  #cat("\n")
  L_m[,i+1] <- this_maf
}
matplot(t(L_m), type = "l")

#It works fine! Ok, so I can now use the influence matrix A to do weighted mafs and drift them
#for ONE maf. Now I just need to figure out how to do this for a list of mafs while simultaniously
#changing pop sizes, then turn that into one function to be able to fully model a single, non-branching
#system. Branches will come later.
#step one is probably to add density dependence, since otherwise the damn thing overflows.











#Density dependence
a <- 2 #set a
R <- 4 #set R
m <- 1 #set m
L <- 10 #set L
n <- 100 #set n
l <- 10 #set number of loci
lmaf <- rep(0.5, length = l) #get random mafs from uniform dist
Tf <- 100

#get the input matrix
L10 <- A(a, R, m, L, n, "both")
L10xs <- unlist(L10[[2]])
L10mat <- unlist(L10[[1]])
#Note: multiplying a matrix with no growth (R = 1) by a different R value is the same as using that
#R value to make the initial matrix. Next question: is getting a pop value by multiplying the resultant
#pop growth by R the same as using the matrix with a higher R? It should be.
L10r1 <- unlist(A(a, 1, m, L, n, "both")[[1]]) #make an R = 1 matrix
N <- ifelse(L10xs < -9 & L10xs > -10, 1, 0) #set N
N2 <- N #set N2
N <- L10mat%*%N #run R = 4 for one gen
N2 <- L10r1%*%N2 #run R = 1 for one gen
N3 <- N2*R #multiply N2 by R = 4...
N3 == N #looks good!

#so multiplying the results by R does work, don't have to include that in the matrix. Can therefore
#still make only one matrix, then weight R by bev-holt.
#bev-holt function to apply to every pop every gen
BH <- function(nt, r, K){
  ntf <- (exp(r)*nt)/(1 + (exp(r) - 1)*(nt/K))
  return(ntf)
}

r = 1.2 #set r
K = 1000 #set K

N <- ifelse(L10xs < -9 & L10xs > -10, 1, 0) #set N
N <- L10r1%*%N #spread
BH(N, r = r, K = K) #grow with bev-holt

#to fully wrap growth in a single loop, then, would look like this...
N <- ifelse(L10xs < -9 & L10xs > -10, 1, 0) #set N
for(i in 1:(Tf-1)){
  N <- BH((L10r1%*%N), r, K)
}
#Great, looks like this works! Density dependence works. Now to try and tie pop growth and Maf changes
#together in one function...










#Doing maf drift and pop growth together:
a <- 2 #set a
R <- 4 #set R
m <- 1 #set m
L <- 10 #set L

#Here was the loop for doing drift...
L_m <- matrix(NA, nrow = length(N), ncol = Tf)
L_m[,1] <- test_maf
for(i in 1:(Tf-1)){
  this_maf <- L_m[,i]
  this_maf <- make_wmafs(in.mafs = this_maf, i.mat = L10mat)
  #cat("Weighted_maf:")
  #print(this_maf)
  for(j in 1:length(this_maf)){#drift in each pop
    this_maf[j] <- make_pop(N[j], 1, this_maf[j])
  }
  #cat("After drift:")
  #print(this_maf)
  #cat("\n")
  L_m[,i+1] <- this_maf
}
matplot(t(L_m), type = "l")

#which I need to combine with the loop for doing pop sizes...
N <- ifelse(L10xs < -9 & L10xs > -10, 1, 0) #set N
for(i in 1:(Tf-1)){
  N <- BH((L10r1%*%N), r, K)
}


#shit, I need to scale maf contributions to both the matrix AND the pop sizes. A spot might have a 
#matrix weight, but if there are no organisms there that shouldn't matter. Likewise, a spot with 
#a huge pop size should matter more. That needs to go into this loop. I should just be able to double
#weight the maf, once by the influence matrix A and once by the pop sizes at all of those points.
#Since the model has spread -> growth, I need to weight by the pop sizes prior to growth


#first, though, I need to modify make_pop to allow it to take a vector of Ns and loci at the same time.


#essentially, I'm dealing with one vector of pop sizes and one for each loci. May as well
#store these as a matrix, but let's first test it with just one loci, storing a matrix for both for
#illustration purposes:

#psuedo for this loop:
#1. Spread the population.
#2. weight starting mafs by influence and by current N.
#3. Grow the pops
#4. use make_pop to drift the input mafs and get output mafs
n <- 100
r <- 1.1
K <- 1000
L10r1 <- unlist(A(a, 1, m, L, n, "both")[[1]]) #make an R = 1 matrix
L10xs <- unlist(A(a, 1, m, L, n, "both")[[2]])

Tf <- 2
N0 <- ifelse(L10xs < -9 & L10xs > -10, 100, 0) #set N0
N <- matrix(NA, nrow = length(N0), ncol = Tf) #initialize pop matrix, rows are site, columns are times
N[,1] <- N0 #add starting pops
L_m <- matrix(NA, nrow = nrow(N), ncol = Tf) #initialize locus matrix, rows are site, columns are time
L_m[,1] <- ifelse(N0 != 0, .45, NA) #add starting mafs, NAs where pop size = 0.


#first need a function to weight by both ns and by influence
testm <- c(rep(.45, 50), rep(.1, 50))
testn <- c(rep(1000, 50), rep(100, 50))

wm.ni <- function(ns, inf, maf){
  w <- inf*ns #weight is the product of the number of individuals and their proximity
  wn <- w/sum(w) #relative weights
  out <- sum(maf*wn) #get the weighted maf for this pop (for the given influence matrix)
  return(out)
}

#...and the function to use this with a whole matrix of weights and return a vector of new mafs.
make.wm.nis <- function(in.mafs, i.mat, ns){
  out <- apply(i.mat, 1, wm.ni, maf = in.mafs, ns = ns)
  return(out)
}

for(i in 1:Tf){
  #spread first:
  N[,i+1] <- L10r1%*%N[,i]
  print(N)
  #weight starting mafs by influence and then current N
  L_m[,i+1] <- make.wm.nis(L_m[,i], L10r1, N[,i+1])
  #grow the pops
  N[,i+1] <- BH(N[,i+1], r, K)
  #drift the mafs
  L_m[,i+1] <- make_pop(n = N[,i+1], l = length(L_m[,i]), maf = L_m[,i])
}



#just realized that this is all probably wrong in the last section. spread should REALLY work by
#determining the number of individuals from each pop that moves using the matrix A, then doing a binom
#test for those. After all, the individuals that move have genotypes drawn randomly from their population.

Tf <- 2
N0 <- ifelse(L10xs < -9 & L10xs > -10, 100, 0) #set N0
N <- matrix(NA, nrow = length(N0), ncol = Tf) #initialize pop matrix, rows are site, columns are times
N[,1] <- N0 #add starting pops
L_m <- matrix(NA, nrow = nrow(N), ncol = Tf) #initialize locus matrix, rows are site, columns are time
L_m[,1] <- ifelse(N0 != 0, .45, NA) #add starting mafs, NAs where pop size = 0.

for(i in 1:(Tf-1)){
  N[,i+1] <- L10r1%*%N[,i]
}

test <- N0*L10r1[1,] #returns the number of individuals which end up in the pop from each other pop
#do random binom tests for these
found.maf <- function(maf, n){
  n <- round(n)
  draw <- numeric(length(maf))
  for (i in 1:length(maf)){
    if(n[i] == 0 | is.na(maf[i]) == TRUE){
      draw[i] <- NA
      next
    }
    #print(n[i])
    #print(maf[i])
    draw[i] <- rbinom(1, 2*n[i], maf[i])
  }
  draw <- sum(draw, na.rm = TRUE)/(sum(n)*2)
  return(draw)
}

found.maf(L_m[,1], test)

#HOWEVER: the weighted maf of the individuals that stay and the individuals which leave 
#one pop should equal the maf of the pop to start with, rather than being simply a binom draw.
#therefore the row of the A matrix matters, and the ns should therefore be calculated in the function
found.maf <- function(maf, i.mat, n0){
  dm <- i.mat*rev(n0) #get the NUMBER of individuals leaving each pop for each other pop
  dm <- dm[dim(dm)[1]:1,] #reverse row order
  dm <- dm[,dim(dm)[1]:1] #reverse column order, dm is now a matrix with the number of inds that
                          #leave each pop for each other pop, with [1,1] the number that stay in 1,
                          #[1,2] the number that leave 1 for 2, and so on.
  #the best thing to do is generate a matrix of mafs for each set of dispersers/remainers using
  #using a binom. Overwrite the stay instances later.
  
  bm <- matrix(NA, dim(dm)[1], dim(dm)[2]) #create an output matrix for mafs after binom
  for(i in 1:nrow(dm)){#loop through rows
    for(j in 1:ncol(dm)){#loop through columns
      cat("i:", i, "j:", j, "dm[i,j]:", dm[i,j], "\n")
      if(dm[i,j] == 0 | is.na(maf[i]) == TRUE){
        bm[i,j] <- NA #if there are no dispersers (in which case the maf[i] should be NA), set the output maf to NA
      }
      else{
        bm[i,j] <- rbinom(1, 2*dm[i,j], maf[i])/2*dm[i,j] #do a binom test to set the maf of the dispersers
      }
    }
  }
  return(bm)
}

maf.o <- .45
n.s <- 10
n.l <- 90
maf.l <- .4

tl <- maf.l*(n.l/(n.l+n.s))
to <- maf.o
ts <- to-tl
maf.s <- ts*(n.l+n.s)/n.s


N0 <- c(500, N0[2:length(N0)])
temp <- L10r1[,1]
test <- L10r1%*%N0
test2 <- L10r1[,1]*rev(N0)
test3 <- L10r1*rev(N0)

round(test[1],7) == round(colSums(test3)[100], 7)
test4 <- test3[dim(test3)[1]:1,]
test5 <- test4[,dim(test3)[1]:1]

N0 <- ifelse(L10xs < -9 & L10xs > -10, 100, 0)
L_m <- ifelse(N0 != 0, .45, NA)
tout <- found.maf(maf = L_m, i.mat = L10r1, n0 = N0)

temp <- tout[1,2:ncol(tout)]
sum((temp*test5[1,2:100])/sum(test5[1,2:100]), na.rm = TRUE)

d <- 2
temp[c(1:(d-1), d+1:length(temp))]





#how do I do binom draws for all individuals at once?
test <- c(10,20,30,10)
nt <- sum(test)
maf <- .45
inds <- rbinom(nt, 2, .45)
omafs <- numeric(length(test))
pos <- 0
for(i in 1:(length(test))){
  omafs[i] <- sum(inds[(pos+1):(pos+test[i])])/(2*test[i])
  pos <- pos + test[i]
}

a <- 2 #set a
R <- 1 #set R
m <- 1 #set m
L <- 10 #set L
n <- 100 #set n
L10 <- A(a, R, m, L, n, "both")
L10xs <- unlist(L10[[2]])
L10mat <- unlist(L10[[1]])
N0 <- ifelse(L10xs < -9 & L10xs > -10, 100, 0) #set N0
L10r1 <- unlist(A(a, 1, m, L, n, "both")[[1]]) #make an R = 1 matrix
L_m <- ifelse(N0 != 0, .45, NA) #add starting mafs, NAs where pop size = 0.
tout <- found.maf(maf = L_m, i.mat = L10r1, n0 = N0)

testL10 <- L10r1
testL10[100,94] <- .000000000001
tout <- found.maf(maf = L_m, i.mat = testL10, n0 = N0)
tout2 <- found.maf(maf = L_m, i.mat = testL10, n0 = N0)
tout[1:5,]
tout2[1:5,]


out <- found.maf(L_m, L10r1, N0)



#ok, so the binom draws aren't working right. Really what I need to do is draw equal to the total population,
#then parse out individuals with allele counts to each other site and stay, since nobody dies here.
#got this working now. found.maf returns the mafs after drift. THIS SHOULD BE ALL I NEED. plug this into a loop
#and I've got non-branching drift for one allele. Do for multiple alleles and I'm golden


#######CURRENT PROBLEMS:
#3.The true pattern should be like this: a) population spreads, doing a binom draw from the maf for each individual.
#b) the population grows with NO NEW BINOM DRAW, since this happens during spread. c) Maf at new location
#calculated, individual data can be dropped. d) the population spreads again.



a <- 2 #set a
R <- 1 #set R
m <- 2 #set m
L <- 10 #set L
n <- 100 #set n
r <- 1.1 #set r
K <- 1000 #set K
Tf <- 100 #set number of gens

L10 <- A(a, R, m, L, n, "both")
L10xs <- unlist(L10[[2]])
L10mat <- unlist(L10[[1]])
N0 <- ifelse(L10xs <= -9 & L10xs >= -10, 100, 0) #set N0
L10r1 <- unlist(A(a, 1, m, L, n, "both")[[1]]) #make an R = 1 matrix

N <- matrix(NA, nrow = length(N0), ncol = Tf) #set N matrix
N[,1] <- N0
L_m <- matrix(NA, nrow = length(N0), ncol = Tf) #set L matrix
L_m[,1] <- ifelse(N0 != 0, .5, NA) #add starting mafs, NAs where pop size = 0.


#This works!
for(i in 1:(Tf-1)){
  L_m[,i+1]<- found.maf(L_m[,i], L10r1, N[,i]) #get new mafs
  N[,i+1]<-L10r1%*%N[,i] #spread pop
  N[,i+1]<-BH(N[,i+1], r, K) #grow pop
}

L_m <- ifelse(is.nan(L_m) == TRUE, NA, L_m)
matplot(t(L_m), type = "l")


#for multiple alleles....
fs <- make_pop(100, l = 10)
l.ms <- matrix(NA, ncol = length(fs), nrow = length(L10xs))
l.ms[1,] <- fs
N0 <- c(100, rep(0, 99))
datm <- cbind(N0, l.ms) #first column is population data, subsequent columns are loci

for(i in 1:(Tf-1)){
  print(i)
  for(j in 2:(ncol(datm))){ #get mafs for each loci
     datm[,j]<- found.maf(datm[,j], L10r1, datm[,1]) #get new mafs
  }
  datm[,1]<-L10r1%*%datm[,1] #spread pop
  datm[,1]<-BH(datm[,1], r, K) #grow pop
}

datm <- cbind(datm, xs = L10xs)

library(ggplot2)
library(reshape2)
colnames(datm) <- c("N0", paste0("loci", 1:10), "xs")
datm_m <- melt(as.data.frame(datm), id.vars = c("N0", "xs"))
fs.df <- data.frame(loci = paste0("loci", 1:10), maf = fs)
colnames(datm_m) <- c("N", "xs", "loci", "maf")
ggplot(datm_m, aes(x = xs, y = maf, color = loci)) + geom_line() + facet_grid(~loci) +
  geom_hline(data = fs.df, aes(yintercept = fs, color = loci))


#---------------------------------------------------------------------------------------------------
#alright, looks good. The model works and is now complete sans-forking, but is very slow. I can make
#this more efficient and wrap it in one function. There is a limit to efficiency, though. Still
#need to do a lot of binom draws.
a <- 2 #set a
m <- 0.5 #set m
L <- 50 #set L
n <- 100 #set n (number of positions along -L to L)
r <- 2 #set r
k <- 3000 #set K
Tf <- 40 #set number of gens
l <- 100 #set number of loci
N0 <- c(100, rep(0, 99)) #set N0

out <- stream.drift(N0, l, a, r, m, L, k, Tf, 0.0002)
library(ggplot2)
library(reshape2)

datm_m <- melt(as.data.frame(out$dat), id.vars = c("N", "xs"))
fs.df <- data.frame(loci = paste0("loci", 1:length(out$imafs)), maf = out$imafs)
colnames(datm_m) <- c("N", "xs", "loci", "maf")

ggplot(datm_m, aes(x = xs, y = maf, color = loci)) + geom_line() + facet_grid(~loci) +
  geom_hline(data = fs.df, aes(yintercept = maf, color = loci))

smaf <- out$imafs
emaf <- out$dat[nrow(out$dat),3:ncol(out$dat)]
comp.df <- cbind(smaf = smaf, emaf = as.numeric(emaf))
comp.df <- as.data.frame(cbind(loci = colnames(out$dat[,3:ncol(out$dat)]), comp.df), stringsAsFactors = FALSE)
comp.df[,2] <- as.numeric(comp.df[,2])
comp.df[,3] <- as.numeric(comp.df[,3])


ggplot(comp.df, aes(x = smaf, y = emaf)) + geom_point() +
  geom_abline(slope = 1, intercept = 0)

temp <- out$dat
for(i in 1:length(smaf)){
  temp[,i+2] <- (out$dat[,i+2]) - smaf[i]
}

temp$mean <- rowSums(temp[,3:ncol(temp)])/(ncol(temp) - 2)

temp_m <- melt(as.data.frame(temp), id.vars = c("N", "xs", "mean"))
colnames(temp_m) <- c("N", "xs", "mean", "loci", "maf")


ggplot(data = temp_m, aes(x = xs, y = maf, fill = loci)) + geom_line(alpha = 0.4) + 
  guides(color = FALSE) + geom_line(aes(x = xs, y = mean), color = "black")

ggplot(temp, aes(x = xs, y = mean)) + geom_point() + geom_hline(yintercept = 0) 

start <- 0
for(i in 1:length(smaf)){
  todo <- (start + 1):(start + (2*L))
  temp_m[todo,6] <- smaf[i]
  start <- todo[length(todo)]
}

ggplot(data = temp_m, aes(x = xs, y = maf, fill = loci, color = smaf)) + geom_line() + 
  xlab("Position") + ylab("Change in Frequency") +
  scale_color_gradient2(low = "blue", mid = "violet", high = "red", midpoint = 0.25)



#----------------------------------------------------------------------------------------------------
a <- 1
m <- 2
L <- 10
n <- 10

#working on forking
A <- function(a, m, L, n){
  K <- function(y, x){#make kernal function
    out <- exp(-abs(y - m - x)/a)/(2*a)
    return(out)
  }
  deltax <- 2*L/n #set range of x
  xs <- seq(-L + deltax/2, L - deltax/2, length = n) #create a list of x values in the range of x
  out <- outer(xs, xs, K)*deltax #create matrix
  return(list(out, xs))
}



K <- function(y, x){#make kernal function, impact of x on y
  out <- exp(-abs(y - m - x)/a)/(2*a)
  return(out)
}

xst <- A(a,m,L,n)[[2]]
test <- A(a,m,L,n)[[1]]
test2 <- K(xst,xst[5]) #impact of 5 on everything
test2 <- round(test2*(2*L/n),3)
plot(y = log(test2), x = xst)
deltax <- (2*L/n)

#this equation below, which is the average of moving x -> y and y -> z produces the movement from
#x -> . 
x <- xst[6]
y <- 0
z <- xst[7]
val <- exp((-abs(y - (m/2) - x)-abs(z - (m/2) - y))/a)/(2*a)
val*deltax
#So... for y -> x then x -> y (say, moving from pos -1, around a point at 0, then back down to -1...)
#this seems to work? Need to make the return matrixes using this kernal instead. ONLY if they 
#go up the other branch (upstream).
#check that they start above, and then move below the point.
#This kernal would be for the ones that then go downstream on the other fork,
#so make sure to halve these influences, since half just continue
#upstream (or whatever the fork weight is). y would be the turnaround point, the position of the fork.
#use this kernal only for the segment of the matrix where, if the river goes A->B OR A->C (fork
#between B and C), B is contributing to C and vice versa.

K2 <- function(x,z){
  out <- exp((-abs(y - (m/2) - x)-abs(z - (m/2) - y))/a)/(2*a)
  return(out)
}
#note that I might have to change the amount that m is divided by by the relative amount transversed
#in each branch. That would make sense. I that case, the kernal would be:
K2 <- function(x,z){
  trav1 <- abs(x - y)
  trav2 <- abs(y - z)
  tot <- trav1 + trav2
  s1 <- trav1/tot
  s2 <- trav2/tot
  out <- exp((-abs(y - (m*s1) - x)-abs(z - (m*s2) - y))/a)/(2*a)
  return(out)
}
#this kernal probably makes more sense, given that it factors in the distance traveled up vs. downstream
#I should probably still check it with someone, though.

#say we have a branch at 0. For anything moving across zero, things are different. If continuing in
#the same direction (up vs. downstream), its easy, simply .5 of the previous influence. If going from
#up to down, its different, need to change the kernal.

example <- matrix(c(1,0.5,0.5,0.5,1,"0.5K2",0.5,"0.5K2",1),nrow = 3, byrow = TRUE)
colnames(example) <- c("A", "B", "C")
rownames(example) <- colnames(example)
example #a matrix with a schematic for a transition matrix with branches C and B, where the number
        #is the weight of the contribution and K2 means to use the alternate Kernal for up->downstream
        #transition. No need to use that with no advection.




#----------------------------------------------------------------------------------------------------
#So, to make this matrix, let's say that we have 200miles of river, with a fork halfway from the start.
#therefore, we'll need to start with the typical matrix. However, we also want to say that there is an
#impermiable barrier upstream at start. So, we'll make a matrix that is much bigger and cleave off the start,
#summing their data into row one. Set L to 150, cleave off the first one hundred

a <- 2 #set a
m <- 0.5 #set m
L <- 150 #set L
n <- 300 #set n (number of positions along -L to L)

A <- function(a, m, L, n){
  K <- function(y, x){#make kernal function
    out <- exp(-abs(y - m - x)/a)/(2*a)
    return(out)
  }
  deltax <- 2*L/n #set range of x
  xs <- seq(-L + deltax/2, L - deltax/2, length = n) #create a list of x values in the range of x
  out <- outer(xs, xs, K)*deltax #create matrix
  return(list(out, xs))
}

mat1 <- A(a, m, L, n)
xs <- mat1[[2]]
mat1 <- mat1[[1]]
#first one hundred columns of this don't matter (contribution of dummy positions to everything else)
mat1 <- mat1[,101:ncol(mat1)]
#sum the first 101 rows, dump the sum into row 101
mat1[101,] <- colSums(mat1[1:101,])
#get rid of first 100 rows
mat1 <- mat1[101:nrow(mat1),]


#Ok, we've dumped the first 100 rows, need to adjust xs...
xs <- xs[101:length(xs)]

#branch is at mile 50. Building the parts of the final matrix. First, from start to mile 50 is the original
#matrix
branch <- 50

#here's the example...
example <- matrix(c(1,0.5,0.5,0.5,1,"0.5K2",0.5,"0.5K2",1),nrow = 3, byrow = TRUE)
colnames(example) <- c("A", "B", "C")
rownames(example) <- colnames(example)
example 

#example["A","A"]
m.a.a <- mat1[1:100, 1:100] #first one-hundred elements, where xs is less than fifty

#example["A","B"]
#needs to be all of the columns where B is contributing to A x0.5
m.a.b <- (mat1[1:100,101:ncol(mat1)])*0.5

#example["A","C"]
#needs to be all of the columns where C (other branch) is contributing to A x0.5. Same as m.a.b if
#branch weight is 50:50
m.a.c <- m.a.b

#example["B","A"]
#all cells where mainstem (A) is contributing to B branch x0.5
m.b.a <- (mat1[101:nrow(mat1),1:100])*.05

#example["B","B"]
#cells where B contributes to B
m.b.b <- mat1[101:nrow(mat1),101:ncol(mat1)]

#example["C", "A"]
# cells where A contributes to C. Same as m.b.a
m.c.a <- m.b.a

#example["C","C"]
#cells where C contributes to C. Same as m.b.b
m.c.c <- m.b.b

#example["B","c"]
#cells where B contributes to C. Needs to use the new kernal! Output is x0.5.
A2 <- function(a, m, L, n, y, xs){
  K2 <- function(x,z){
    trav1 <- abs(x - y)
    trav2 <- abs(y - z)
    tot <- trav1 + trav2
    s1 <- trav1/tot
    s2 <- trav2/tot
    out <- exp((-abs(y - (m*s1) - x)-abs(z - (m*s2) - y))/a)/(2*a)
    return(out)
  }
  deltax <- 2*L/n #set range of x
  out <- outer(xs, xs, K2)*deltax #create matrix
  return(list(out, xs))
}

b.xs <- xs[101:length(xs)]
m.b.c <- A2(a,m,L,n, y = branch, xs = b.xs)*0.5

#example["C","B"]
#same as m.b.c
m.c.b <- m.b.c

#That's all nine. Need to put them together.
mat2 <- rbind(cbind(m.a.a, m.a.b, m.a.c), cbind(m.b.a, m.b.b, m.b.c), cbind(m.c.a, m.c.b, m.c.c))

#ok, let's try stream.drift on this.
n0 <- c(100, rep(0,299)) #set n0
r <- 2 #set r
k <- 3000 #set K
Tf <- 40 #set number of gens
l <- 10 #set number of loci
xs2 <- c(xs[1:100], xs[101:length(xs)], xs[101:length(xs)]) #set xs

out <- stream.drift(n0, l, NULL, r, NULL, NULL, k, Tf, in.mat = mat2, in.xs = xs2)
out$dat$branch <- c(rep("mainstem", 100), rep("A", 100), rep("B", 100))




#and plot some stuff......
library(ggplot2)
library(reshape2)
smaf <- out$imafs


#plots of final change in allele frequencies
temp <- out$dat
for(i in 1:length(smaf)){
  temp[,i+2] <- (out$dat[,i+2]) - smaf[i]
}
temp_m <- melt(as.data.frame(temp), id.vars = c("N", "xs","branch"))
colnames(temp_m) <- c("N", "xs", "branch", "loci", "maf")

start <- 0
for(i in 1:length(smaf)){
  todo <- (start + 1):(start + (2*L))
  cat(todo[1], todo[length(todo)], "\n")
  temp_m[todo,6] <- smaf[i]
  cat(temp_m[todo[1],6], temp_m[todo[length(todo)],6], "\n")
  start <- todo[length(todo)]
}
colnames(temp_m) <- c("N", "xs", "branch", "loci", "maf", "smaf")

ggplot(data = temp_m, aes(x = xs, y = maf, color = loci, linetype = branch)) + geom_line() + 
  xlab("Position") + ylab("Change in Frequency") +
  geom_hline(yintercept = 0)

ggplot(data = temp_m, aes(x = xs, y = maf, fill = loci, color = smaf, linetype = branch)) + 
  geom_line() + 
  scale_color_gradient2(low = "blue", mid = "violet", high = "red", midpoint = 0.25) +
  xlab("Position") + ylab("Change in Frequency") +
  geom_hline(yintercept = 0)


#graph of final allele frequencies
out_m <- melt(as.data.frame(out$dat), id.vars = c("N", "xs","branch"))
colnames(out_m) <- c("N", "xs", "branch", "loci", "maf")

start <- 0
for(i in 1:length(smaf)){
  todo <- (start + 1):(start + (2*L))
  cat(todo[1], todo[length(todo)], "\n")
  out_m[todo,6] <- smaf[i]
  cat(out_m[todo[1],6], out_m[todo[length(todo)],6], "\n")
  start <- todo[length(todo)]
}
colnames(out_m) <- c("N", "xs", "branch", "loci", "maf", "smaf")

ggplot(data = out_m, aes(x = xs, y = maf, color = loci, linetype = branch)) + geom_line() + 
  xlab("Position") + ylab("Allele Frequency") +
  geom_segment(data = out_m, aes(y =  smaf, yend = smaf, color = loci), x = -75, xend = -50) +
  geom_vline(xintercept = 50)

ggplot(data = out_m, aes(x = xs, y = maf, fill = loci, color = smaf, linetype = branch)) + 
  geom_line() + 
  scale_color_gradient2(low = "blue", mid = "violet", high = "red", midpoint = 0.25) +
  xlab("Position") + ylab("Ending in Frequency") +
  geom_hline(yintercept = 0)

#plot endpoint drift vs starting in point
ends <- out_m[out_m$xs == max(out_m$xs) | out_m$xs == max(out_m[out_m$branch == "mainstem",]$xs),]
ggplot(ends, aes(x = smaf, y = maf, shape = branch, color = branch)) + geom_point(size = 3) + 
  geom_abline(intercept = 0, slope = 1) + xlab("Starting Frequency") + ylab("Ending Frequency")

ends$branch <- as.factor(ends$branch)
mod.A <- with(ends, lm(maf~smaf + branch))
mod.B <- with(ends, lm(maf~smaf + relevel(branch, ref = "B")))
mod.m <- with(ends, lm(maf~smaf + relevel(branch, ref = "mainstem")))
mod.simp <- with(ends, lm(maf~smaf))
mod.change <- with(ends, lm(abs((maf-smaf))~smaf))

#try just doing growth
N <- matrix(NA, 300, Tf)
N[,1] <- n0

for(i in 1:(Tf-1)){
  N[,i+1] <- BH(mat2%*%N[,i], r, k)
}
N <- as.data.frame(N)
colnames(N) <- 1:Tf
N$pos <- temp$xs
N$branch <- temp$branch

N_m <- melt(N, id.vars = c("pos", "branch"))
colnames(N_m) <- c("pos", "branch", "yr", "N")



g_legend<-function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)}

Npm <- ggplot(N_m[N_m$branch == "mainstem",], aes(x = pos, y = N, color = yr)) + geom_line() +
  scale_color_discrete(name = "Year")
legend <- g_legend(Npm)
Npm <- ggplot(N_m[N_m$branch == "mainstem",], aes(x = pos, y = N, color = yr)) + geom_line() + theme(legend.position="none") + xlab("Position")
NpA <- ggplot(N_m[N_m$branch == "A",], aes(x = pos, y = N, color = yr)) + geom_line() + theme(legend.position="none") + xlab("Position") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
NpB <- ggplot(N_m[N_m$branch == "B",], aes(x = pos, y = N, color = yr)) + geom_line() + theme(legend.position="none") + xlab("Position")



library(gridExtra)
Npall <- grid.arrange(Npm, NpA, NpB, legend, ncol = 7, layout_matrix = cbind(c(1,1),c(1,1), c(1,1), c(2,3), c(2,3), c(2,3), c(4,4)))
#purty!













####################################################################################################################
#working on function to create matrix from a df of branches

#a): <branch.id><start.pos><start.node><start.weight><end.pos><end.node><end.weight>


library("igraph")
set.seed(1702)


bs <- data.frame(branch.id = c("001", "002", "003", "004", "005", "006"), 
                 length = c(100, 25, 100, 75, 25, 20),
                 start.node = c("C", "B", "A", "D", "E", "D"),
                 start.weight = c(1, .5, 1, 1, 1, 1),
                 end.node = c("A", "A", "F", "B", "B", "C"),
                 end.weight = c(1,1,.5,.2,.8,1), stringsAsFactors = FALSE)




bsm <- as.matrix(bs[,c(3,5)])
rownames(bsm) <- bs[,1]

nbel <- graph_from_edgelist(bsm, directed = FALSE)
all_simple_paths(nbel, "D", "B")
E(nbel)$name <- bs[,1]
E(nbel)$s.weight <- bs[,4]
E(nbel)[E(nbel)$name == 1]$s.weight

temp <- nbel + vertex(c("s", "e")) + 
  edges(c("D", "s"), c("B", "s"), c("A", "e"), c("F", "e")) -
  edges(get.edge.ids(nbel, c("D", "B")), get.edge.ids(nbel, c("A", "F")))
paths <- all_simple_paths(temp, "s", "e")
verts1 <- paths[[1]][-c(1, length(paths[[1]]))]
verts2 <- paths[[2]][-c(1, length(paths[[2]]))]

unique(c(bs$start.pos[bs$start.node %in% V(nbel)$name[verts2][2]], 
         bs$end.pos[bs$end.node %in% V(nbel)$name[verts2][2]]))


inlist <- c()
for(i in 1:length(verts2)){
  inlist <- c(inlist, unique(c(bs$start.pos[bs$start.node %in% V(nbel)$name[verts2][i]], 
    bs$end.pos[bs$end.node %in% V(nbel)$name[verts2][i]])))
}

inlist <- c()
for(i in 1:(length(verts2)-1)){
  inlist <- c(inlist, 
              get.edge.ids(nbel, c(V(nbel)$name[verts2][i], V(nbel)$name[verts2][i+1])))
}

get.edge.ids(nbel, c(V(nbel)$name[verts2][-length(verts2)], V(nbel)$name[verts2][-1]))

test1 <- c(1,2,3,4)
test2 <- c(5,6,7,8)
outer(test1, test2) #in outer product, rows are the first entry, columns are the second
#since we are going FROM sport marked by row TO spot marked by column, this would be:
# outer(xs.s, xs.e)
a <- .5
m <- 1
L <- 10
n <- 5
#out <- A(a, m, L, n, "both")

#K(-8, -4)*2*L/n
#confirming, yes, looks correct

#continuing with going from branch 4 to branch 2, passing through 6 and 1
dx <- 1
xs.s <- seq(-bs[4,2]/2 + dx/2, bs[4,2]/2 - dx/2, by = dx) #calculate the xs for this branch (ending)
xs.e <- seq(dx/2, bs[3,2] - dx/2, by = dx) #ending xs, most upstream at 0



A3 <- function(a,m,dx,y,xs.s,xs.e){
  K3 <- function(x,z){
    y <- c(x,y,z)
    part <- numeric(length(y) - 1)
    cat("This y:", y, "\n")
    for(i in 1:(length(y) - 1)){
      #print(i)
      trav <- abs(y[i] - y[i+1])
      s <- trav/sum(abs(y[-1] - y[-length(y)]))
      #print(sum(abs(y[-1] - y[-length(y)])))
      #cat("y[i], y [i+1]:", y[i], y[i+1], "\n")
      part[i] <- -abs(y[i+1] - (m*s) - y[i])
      #print(part)
    }
    out <- exp(sum(part)/a)/(2*a)
    print(out)
    return(out)
  }
  out <- outer(xs.s, xs.e, Vectorize(K3))*dx #create matrix
  return(out)
}

out <- A3(a, m, dx, y, xs.s[1:2], xs.e[1:2])
K3(xs.s[1], xs.e[1])


test1 <- c(1,2,3)
test2 <- c(5,6)
outer(test1, test2, "+")

y <- ifelse(bs[4,3] == V(nbel)$name[verts2[1]], 0, bs[4,2])
w <- numeric(length(verts2))
lvert <- "B"
for (i in 1:length(inlist)){
  tvert <- V(nbel)$name[verts2[i]]
  y[i+1] <- ifelse(bs[inlist[i],3] == tvert, y[i] + bs[inlist[i],2], y[i] - bs[inlist[i],2])
  wi <- ifelse(bs[inlist[i],3] == tvert, bs[inlist[i],4], bs[inlist[i],6])
  forks <- bs[which((bs[,3] == tvert | bs[,5] == tvert) & 
                      (bs[,3] != lvert & bs[,5] != lvert)),]
  #print(forks)
  w[i] <- wi/sum(ifelse(forks[,3] == tvert, forks[,4], forks[,6]))
  lvert <- tvert
}
w[length(w)] <- ifelse(bs[3,3] == V(nbel)$name[verts2[length(verts2)]],bs[3,4], bs[3,6])
w[length(w)] <- w[length(w)]/2 #shorthand, process to do this in function

#giving vertex its weight attributes, named by destination node
#for (h in 1:nrow(bs)){
#  set_vertex_attr(nbel,"wt", V(nbel)[V(nbel)$name == bs[h,3]], list(c(paste0("wt",bs[h,5]) = bs[h,4])))
#}

#paste0("wt",bs[h,5])
#set_vertex_attr(nbel,paste0("wt.",bs[1,5]), V(nbel)[V(nbel)$name == bs[1,3]],bs[1,4])
#bs[bs[,3] == V(nbel)$name[V(nbel)[nbel[1] == 1][1]],5]
#V(nbel)$name[V(nbel)[nbel[1] == 1][1]]

temp <- bs[which((bs[,3] == "A" | bs[,5] == "A") & (bs[,3] != "C" & bs[,5] != "C")),]

ws <- ifelse(temp[,3] == "A", temp[,4], temp[,6])

test1 <- matrix(c(2,3,4,5,6,7), nrow = 2)
test2 <- matrix(c(4,5,6,7,8,9), nrow = 2)
testl <- list(test1, test2)


#test it
test_out <- create.matrix(bs, a, 5, dx = 1, F)
