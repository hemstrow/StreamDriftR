k
y%*%x

L<-100 #100 meters
n<- 30 #number of spots
m<- 1 #mean location they disperse to
a<- 5 #width of the dispersal 
K <- function(y, x){#make kernal function
  out <- exp(-abs(y - m - x)/a)/(2*a)
  return(out)
}
deltax <- 2*L/n #set range of x
xs <- seq(-L + deltax/2, L - deltax/2, length = n) #create a list of x values in the range of x
out <- outer(xs, xs, K)*deltax #create matrix

#A function to calculate population growth based on the beverton-holt model:
#Inputs: nt: starting population
#        r: intrinsic growth rate
#        K: carrying capacity

BH <- function(nt, r, K){
  ntf <- (exp(r)*nt)/(1 + (exp(r) - 1)*(nt/K))
  return(ntf)
}
#spread and growth in one generation
n0<- numeric(n) #n at time 0
n0[1]<-100
spread<-n0%*%out#spread is number of individuals in each location
growth<- BH(spread, 1,3000)
r<-1 #growth rate
K<-3000 #carrying capacity
##BH model 10 generations
N<-numeric(10)
N[1]<-100
for(i in 2:length(N)){
  N[i]<-BH(N[i-1],r,K)
  deltax <- 2*L/n #set range of x
  xs <- seq(-L + deltax/2, L - deltax/2, length = N) #create a list of x values in the range of x
  out <- outer(xs, xs, K)*deltax
##*deltax ##individuals can disperse anywhere along the 200 meter stretch
}


