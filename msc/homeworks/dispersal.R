# Input: arugments for A, N vector, and arguments for BH
# output: vector of fish numbers after movement
out <- A(1, -5, 50, 99)
N <- numeric(99)
N[1] <- 100
BH(N%*%out, 3, 100)

# Next step: Have this run for g generations, and save the output from each generation in a matrix! (loops!)
disperse_and_grow <- function(N, m, n, L, K, r, g){
  
}



# loops
x <- 101:200
for(i in 2:100){
  x[i-1] <- x[i-1] + x[i]
}
