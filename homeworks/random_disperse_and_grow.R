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
  output[1,] <- BH(colSums(mc2d::rmultinomial(length(N), N, prob = out)), r, K)
  
  for(i in 2:g){
    output[i,] <- colSums(mc2d::rmultinomial(length(N), output[i-1,], prob = out))
    output[i,] <- BH(output[i,], r, K)
  }
  return(output)
}

matplot(out)
