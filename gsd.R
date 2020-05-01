#afs x position
#' Spread individuals along a one-dimensional system, tracking genetic drift
#' 
#' @param L numeric, one half the available range (modeled from -L to L)
#' @param n numeric, the number of populations within range
#' @param m numeric, average displacement distance
#' @param a numeric, 2a^2 is the variance in displacement
#' @param N numeric vector of length n, number individuals at each position 
#' @param afs numeric matrix, row for each loci, column for each position
#' @param loci numeric, default. Number of loci
#' @param r numeric, growth rate
#' @param K numeric, carrying capacity
#' @param gens numeric vector, number of generations to loop through


gsd<- function(L, n, m, a, N, afs, loci, r, K, gens){
  
  #1) spread
  #N%*%dmat
  #rows = source position
  #cols = destination position
  #ex: pos[3,4] means # fish going from 3 to 4
  
  dmat<- A(L/2, n, m, a)
  spread<- function(N, dmat){
    #error: pos only has 1 row--> why?
    #do this as multinomial instead--> rmultinom
    #we will treat dmat as the prob a fish goes to location, multinom tells how many fish went
    #vectorize this similar to drift function
    pos<- mc2d::rmultinomial(nrow(dmat[[1]]), N, dmat[[1]])
    #pos<- rmultinom(N, nrow(dmat[[1]]), dmat[[1]])
    #pos<- N*t((dmat)[[1]])
    return(pos)
  }
  
  #2) find afs
  #2ii) draw alleles at each loci
  #rbinom: 
  #obs = #sites * #loci
  #  sites = #r * #c of pos
  #size = pos, repeated 2 times for each loci, in the order where we have loci reps for each value
  #  if 100 fish stay in site 1 and theres 10 loci, the first 10 values should be 100
  #  if loci = 10, the pos should be repeated 10 times
  #prob = prob of getting minor allele for groups of fish at each loci from source pop
  
  sites<- (nrow(dmat[[1]]))*ncol(dmat[[1]])

  drift<- function(sites, pos, afs, loci){
    
    #repeat afs for each source location, for each loci
    
    repafs<- afs[,rep(1:nrow(pos), each = loci)]
    
    #number of minor alleles
    num_minor<- rbinom(sites*loci, round(rep(pos, each = loci)), repafs)
    #if(any(is.na(num_minor))){
     # browser()
   # }
    
    #data frame of dest #, loci #, and num of minor alleles
    mindf<- data.frame(dest = rep(rep(1:ncol(pos), each = loci), times = ncol(pos)), loci = rep(1:loci, times = ncol(pos)^2), num_minor = num_minor)
    
    #num allele at each dest at each loci
    numminloci<- reshape2::dcast(mindf, dest~loci, fun.aggregate = sum)
    
    #afs at each dest at each loci
    num_fish<- colSums(pos)
    dest_afs<- t(numminloci[,-1]/(2*num_fish))
    return(dest_afs)
  }
  
  #grow
  BH <- function(nt, r, K){
    ntf <- (exp(r)*nt)/(1 + (exp(r) - 1)*(nt/K))
    return(ntf)
    
  }
  
  
  #3) loop through multiple gens
  #3i) initialize array with dest afs for each gen
  
  result<- array(NA, dim = c(loci,ncol(dmat[[1]]), gens))
  dest_afs<- afs
  #3ii) do afs for each gen
  
  for(i in 1:gens){
    
    
    
    pos<- spread(N, dmat)
    
    #do afs
    dest_afs<- drift(sites, pos, dest_afs, loci)
    
    #grow
    N<- BH(colSums(pos), r, K)
    #browser()
    
    #put afs into array for gen
    
    result[,,i]<- dest_afs
    
  }
  
  return(result)
}
