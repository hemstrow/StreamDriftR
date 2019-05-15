
the.function <- function(L, n, m, a, r, K, gen){
  browser()
  #A function to calculate population growth based on the beverton-holt model:
  #Inputs: nt: starting population
  #        r: intrinsic growth rate
  #        K: carrying capacity
  A <- function(L, n, m, a){
    browser()
    K <- function(y, x){#make kernal function
      out <- exp(-abs(y - m - x)/a)/(2*a)
      return(out)
    }
    deltax <- 2*L/n #set range of x
    xs <- seq(-L + deltax/2, L - deltax/2, length = n) #create a list of x values in the range of x
    out <- outer(xs, xs, K)*deltax #create matrix
    return(out)
  }
  
  BH <- function(nt, r, K){
    ntf <- (exp(r)*nt)/(1 + (exp(r) - 1)*(nt/K))
    return(ntf)
  }
  
  out <- A(L, n, m, a)
  
  
  N<-numeric(n)
  N[1]<-100
  
  m <- matrix(NA, gen, n)
  m[1,] <- N
  
  
  for(i in 2:gen){
    m[i,] <- m[i - 1,]%*%out#spread is number of individuals in each location
    m[i,] <- BH(m[i,], r, K)
  }
  return(m)
}
