# Input: arugments for A, N vector, and arguments for BH
# output: vector of fish numbers after movement
out <- A(1, -5, 50, 99)
N <- numeric(99)
N[1] <- 100
BH(N%*%out, 3, 100)

# Next step: Have this run for g generations, and save the output from each generation in a matrix! (loops!)

#' Disperse and grow populations
#' 
#' Models the dispersal and growth of populations along
#' a one-dimensional surface at discritized postions.
#' Growth is modeled usign the Beverton-Hold (BH) growth
#' model. Disperal is modeled using a laplacian kernal.
#' 
#' @param a numeric, variance in dispersal
#' 
#' @return Matrix containing the number of individuals at each position (colunms) and
#' generation (rows).
#' 
#' @example
#' out <- disperse_and_grow(1, c(100, rep(0, 98)), -5,  50, 3000, 3, 100)
#' matplot(out)
disperse_and_grow <- function(a, N, m, L, K, r, g){
  A <- function(a, m, L, n, output = "matrix"){
    K <- function(y, x){#make kernal function
      out <- exp(-abs(y - m - x)/a)/(2*a)
      return(out)
    }
    deltax <- 2*L/n #set range of x
    xs <- seq(-L + deltax/2, L - deltax/2, length = n) #create a list of x values in the range of x
    
    out <- outer(xs, xs, K)*deltax #create matrix
    if(output == "matrix"){
      return(out)
    }
    else if (output == "xs"){
      return(xs)
    }
    else if (output == "both"){
      return(list(out, xs))
    }
    else{
      warning("Please specify an output format.")
      stop
    }
  }
  BH <- function(nt, r, K){
    ntf <- (exp(r)*nt)/(1 + (exp(r) - 1)*(nt/K))
    return(ntf)
  }
  
  n <- length(N)
  out <- A(a, m, L, n)

  output <- matrix(NA, nrow = g, ncol = n)
  output[1,] <- BH(N%*%out, r, K)
  
  for(i in 2:g){
    output[i,] <- BH(output[i-1,]%*%out, r, K)
  }
  return(output)
}


out <- disperse_and_grow(1, c(100, rep(0, 98)), -5,  50, 3000, 3, 100)


# loops
x <- 101:200
for(i in 2:100){
  x[i-1] <- x[i-1] + x[i]
}

