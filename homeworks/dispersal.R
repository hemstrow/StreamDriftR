# Input: arugments for A, N vector, and arguments for BH
# output: vector of fish numbers after movement
out <- A(1, 5, 50, 98)
N <- numeric(98)
N[1] <- 100
BH(N%*%out, 3, 100)

# Next step: Have this run for g generations, and save the output from each generation in a matrix! (loops!)
