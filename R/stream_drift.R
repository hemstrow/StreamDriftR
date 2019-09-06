#funciton to make simulated snp data for a population of a specificed n (number ind), l (number loci),
#maf (minor allele frequency), and mafsd (standard deviation of minor allele frequency).
#If mafsd = FALSE, then maf should be a vector of obseved mafs  calculated from data.
#If maf = NA, a uniform distribution of mafs between maf.min and .5 is used.
#If maf is a numeric vector, it uses that vector.
#Other arguments:
#maf.min: The minimum minor allele frequency to use if maf is not NA.
#pop.out: If FALSE, function outputs only the mafs for each loci after resampling. If true,
#         outputs a simulated population.
#mu:      Mutation rate, minor to major and major to minor only.
#Currently assumes no LD. Note, if mafsd is not uniform and maf is a number (mafs generated
#from a normal distribution), mean mafs will be biased high (away from 0), since mafs of less than 0
#are rejected. This bias is higher with higher sds.
#Current: if n is a vector of length != 1, uses each n with each corresponding maf. This is essentially
#for one locus with multiple mafs in multiple pops, needs to be looped for more loci

make_pop <- function(n, l, maf = NA, mafsd = "uniform", maf.min = 0.05, mu = 0, pop.out = FALSE){
  n <- round(n) #round n, since you can't take a binomial of a decimal
  #print(n)
  #print(maf)
  ##make vector of mafs for each loci based on truncated normal distribution
  if(is.numeric(maf) == TRUE){
    if(length(maf) != l){
      warning("Length of maf vector not equal to the number of loci, sampling from maf vector")
      lmaf <- sample(maf, l, replace = TRUE)
    }
    else{
      lmaf <- maf #set mafs equal to the input vector if given
      #cat("set lmaf to maf:", lmaf, "\n")
    }
  }
  else if (mafsd == "uniform"){
    #cat("Using uniform maf.\n")
    lmaf <- runif(l, maf.min, 0.5) #gets random lmafs between maf.min and 0.5
    #print(lmaf)
  }
  else if (is.numeric(mafsd) == TRUE){
    #print(maf)
    if(maf <= 0 || maf >= 1){
      stop("maf is not between 0 and 1!")
    }
    lmaf <- rnorm(n = l, mean = maf, sd = mafsd) #generate mafs from mafsd and mean maf
    lmaf <- lmaf[lmaf > 0 & lmaf < 1] #truncate
    while(length(lmaf) < l){ #while lmaf stays too short...
      lmaf <- c(lmaf, rnorm(n = l, mean = maf, sd = mafsd)) #add more draws
      lmaf <- lmaf[lmaf > 0 & lmaf < 1] #truncate
    }
    lmaf <- lmaf[1:l] #shorten it down to length l
  }
  else{
    warning("Please specify mafsd")
    stop()
  }
  
  ##use vector of mafs to run the binomial function and generate either mafs for the next generation
  ##or a simulated population.
  #write function to generate n genotypes according to data, if pop.out = TRUE
  if(pop.out == TRUE){
    gen_genos <- function(tmaf){
      #print(tmaf)
      draws <- rbinom(n, 2, tmaf)
      draws_geno <- ifelse(draws == 0,"0101", ifelse(draws == 1, "0102", "0202"))
      return(draws_geno)
    }
    
    #apply genotypes function
    gens <- lapply(FUN = gen_genos, lmaf) #create all data
    gens <- as.data.frame(matrix(unlist(gens), nrow = l, ncol = n, byrow = TRUE), stringsAsFactors = FALSE) #put data in data.frame
    colnames(gens) <- unlist(lapply("ind", paste0, 1:n))
    rownames(gens) <- unlist(lapply("loci", paste0, 1:l))
  
    #output data
    return(gens)
  }
  
  #do random draws, create vector of mafs if pop.out = FALSE
  else{
    if(length(n) == 1){ #if only one pop size was passed...
      #function which does 2*n draws with prob tmaf. Successes are instances where q is passed.
      gen_nmaf <- function(tmaf){
        #print(tmaf)
        draw <- rbinom(1, 2*n, tmaf)/(n*2)
        return(draw)
      }
      
      #get new mafs for each loci, output
      l.n.maf <- unlist(lapply(FUN = gen_nmaf, lmaf))
      return(l.n.maf)
    }
    else{
      gen_nmaf <- function(maf, n){
        draw <- numeric(length(maf))
        for (i in 1:length(maf)){
          draw[i] <- rbinom(1, 2*n[i], maf[i])/(n[i]*2)
        }
        return(draw)
      }
      l.n.maf <- gen_nmaf(maf,n)
      return(l.n.maf)
    }
  }
}

#Function to create matrix of influence of inidividuals in one location to the next.
#parameters:
# a: 2a^2 is the variance in displacement
# R: growth rate
# m: average displacement distance
# L: one half the available range (modeled from -L to L)
# n: the number of populations within range -L:L to examine
# out: if "matrix," returns the matrix, if "xs", returns the series of x values for which the matrix
#      is constructed, if "both," does a list where the matrix is element 1 and xs is element 2.
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


#Function to get a single maf value (output), determined a vector of mafs (mafs) and a vector of
#their weights (weights). Weight in this case is highest for the input maf
#from the same pop that is returned. Automatically calculates relative weights.
w.maf <- function(weights, mafs){
  out <- sum(mafs*weights/sum(weights))
  return(out)
}

#Function which applies w.maf to a matrix. Takes, as in w.maf, a vector of mafs (in.mafs),
#and a matrix (i.mat) where each row is a set of weights. Returns a vector of weighted mafs.
#Assuming the input matrix is a position influence vector (as in a leslie matrix, or one produced
#by the A function), this in essence determines the weighted input mafs for a single loci
#for generation t + 1 given those mafs in generation t. These output mafs still need to be allowed to
#drift using make_pop
make_wmafs <- function(in.mafs, i.mat){
  out <- apply(i.mat, 1, w.maf, mafs = in.mafs)
}


#function which weights by both pop size and by influence of nearby mafs.
#inputs: ns: vector of population sizes
#        inf: matrix of influences based on kernal e.g. laplacian (A)
#        maf: vector of mafs for each pop, in the same order as the pops in inf and ns.
wm.ni <- function(ns, inf, maf){
  w <- inf*ns #weight is the product of the number of individuals and their proximity
  wn <- w/sum(w) #relative weights
  out <- sum(maf*wn) #get the weighted maf for this pop (for the given influence matrix)
  return(out)
}


#Function which applies in.mafs across a matrix of different influence weights. The returned
#vector is the maf for each position, weighted by the pop size*proximity influence of all mafs for the
#loci
make.wm.nis <- function(in.mafs, i.mat, ns){
  out <- apply(i.mat, 1, wm.ni, maf = in.mafs, ns = ns)
  return(out)
}


#A function to calculate population growth based on the beverton-holt model:
#Inputs: nt: starting population
#        r: intrinsic growth rate
#        K: carrying capacity
BH <- function(nt, r, K){
  ntf <- (exp(r)*nt)/(1 + (exp(r) - 1)*(nt/K))
  return(ntf)
}



#function which calculates mafs after dispersal by calculating the number dispersing from a transition
#matrix A (i.mat), starting mafs (maf), and starting pop sizes (n0). The dispersers are caclulated
#as in %*%, but without summing down columns. Mafs are calculated by doing binom draws for each individual,
#dispersing the required number into each spot, then calculating the maf of the whole column (pop),
#outputs mafs.
found.maf <- function(maf, i.mat, n0){
  dm <- t(i.mat)*n0 #get the NUMBER of individuals leaving each pop for each other pop
                    #with [1,1] the number that stay in 1,
                    #[1,2] the number that leave 1 for 2, and so on.
  dm <- round(dm) #round the pops
  bm <- matrix(NA, dim(dm)[1], dim(dm)[2]) #create an output matrix for mafs after binom
  
  for(i in 1:nrow(dm)){#loop through rows
    nds <- dm[i,]#nds is a vector number of pops who went to each spot
    nt <- sum(nds) #the total number of individuals
    if(nt == 0){ #no binoms to do if there are no pops that disperesed, just return NA for the row
      bm[i,] <- NA
    }
    else if(is.nan(maf[i]) == TRUE){
      #cat("we have a problem.\n")
      bm[i,] <- NA #in instances where individuals are dispersing from a pop that wasn't assigned
                   #a maf, this implies that the number of individuals who went into the location
                   #last gen was less than 1 (rounded out). Since you can't have 0 individuals move,
                   #their mafs were zero, and nothing can actually move out this gen.
                   #Therefore, the maf contribution of these individuals is NA, and should be ignored
                   #when calculaing weighted mafs.
      dm[i,] <- 0 #likewise, reset the number of inds moving to zero to avoid messing up weights later
      #browser()
    }
    else{
      inds <- rbinom(nt, 2, maf[i]) #do a draws for each individual from the input maf
      pos <- 0 #set starting position in the inds vector to zero
      for(j in 1:(length(nds))){ #partition individuals form the inds vector into each location
        if(nds[j] == 0){bm[i,j] <- NA} #if no individuals disperesed here, set the allele count as NA and move on
        else{
          bm[i,j] <- sum(inds[(pos+1):(pos+nds[j])]) #get the allele count for each spot
          pos <- pos + nds[j]
        }
      }
    }
  }
  mafs <- colSums(bm, na.rm = TRUE)/(2*colSums(dm, na.rm = TRUE))
  #print(mafs)
  return(mafs)
}


#The full function to simulate drift along a single unbranching stream. Returns a list where the
#first element $dat is the output data, with final N in the first column, the position (xs) in the
#second, and the final mafs in the rest, and the second element $imafs is the initial minor allele
#frequency prior to drift.
#inputs: 
# n0: initial populations
# l: number of loci
# a: 2a^2 is the variance in displacement for dispersal
# r: intrinsic growth rate
# m: displacement distance
# L: one half the available range (modeled from -L to L)
# k: carying capacity for each location. If a single value, k is the same for each location. If a
#    vector, each entry is the k for each location, in the same order as the matrix
# Tf: time to run drift simulation for
# maf.min: minimum allele frequency for randomly generated alleles.
# maf.max: maximum allele frequency for randomly generated alleles.
# in.mat: Optional input matrix if one is not to be created. In this case, a, m, and L are not used
#         and should probably be set to NULL.
# in.xs: An input set of positions to use. Needed if in.mat is defined. Just needed for 
#        printing to xs column in output.
# in.freqs: Optional input set of allele frequencies to be used instead of randomly drawn ones.
stream.drift <- function(n0, l, a = NULL, r, m = NULL, L = NULL, k, Tf, maf.min = 0.05, maf.max = 0.5,
                         in.mat = NULL, in.xs = NULL, in.freqs = NULL){
  #A function to calculate population growth based on the beverton-holt model:
  #Inputs: nt: starting population
  #        r: intrinsic growth rate
  #        K: carrying capacity
  BH <- function(nt, r, K){
    ntf <- (exp(r)*nt)/(1 + (exp(r) - 1)*(nt/K))
    return(ntf)
  }
  
  
  #if no in.mat supplied
  if(is.null(in.mat) == TRUE){
    #a function:
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
    #get matrix
    i.mat <- unlist(A(a, m, L, length(n0))[[1]])
    xs <- unlist(A(a, m, L, length(n0))[[2]])
  }
  
  #if a matrix is supplied:
  else{
    i.mat <- in.mat
    xs <- in.xs
    if(is.null(xs) == TRUE){
      warning("An input xs must be provided if an input i.mat is also given.")
      stop()
    }
  }
  #get and set initial mafs
  if(is.null(in.freqs) == TRUE){
    imaf <- runif(l, maf.min, maf.max) #generate initial mafs
  }
  else{
    if(l != length(in.freqs) | is.numeric(in.freqs) == FALSE){stop("l must equal length of in freqs, which must be numeric!")}
    imaf <- in.freqs
  }
  
  #get and set initial  n0.
  datm <- matrix(NA, length(xs), (2+l)) #initialize output matrix
  datm[,1] <- n0 #set n0 in output matrix
  datm[,2] <- xs #set xs in output matrix
  
  # add allele frequencies where n0 is not zero
  datm[which(n0 != 0),-c(1:2)] <- imaf
  
  for(t in 1:Tf){ #for each time point
    if(t %% 10 == 0){cat("Time:", t, "\n")}
    #set up dm matrix:
    dm <- t(i.mat)*datm[,1] #get the NUMBER of individuals leaving each pop for each other pop
    #with [1,1] the number that stay in 1,
    #[1,2] the number that leave 1 for 2, and so on.
    dm <- round(dm) #round the pops
    
    
    #do all of the binoms for each locus
    for(h in 1:l){ #for each loci...
      bm <- matrix(NA, dim(dm)[1], dim(dm)[2]) #create an output matrix for mafs after binom
      maf <- datm[,(h+2)] #set mafs to work with
      
      for(i in 1:nrow(dm)){#loop through rows of dm
        nds <- dm[i,]#nds is a vector number of pops who went to each spot
        nt <- sum(nds) #the total number of individuals
        if(nt == 0){ #no binoms to do if there are no pops that disperesed, just return NA for the row
          bm[i,] <- NA
        }
        else if(is.nan(maf[i]) == TRUE){
          #cat("we have a problem.\n")
          bm[i,] <- NA #in instances where individuals are dispersing from a pop that wasn't assigned
          #a maf, this implies that the number of individuals who went into the location
          #last gen was less than 1 (rounded out). Since you can't have 0 individuals move,
          #their mafs were zero, and nothing can actually move out this gen.
          #Therefore, the maf contribution of these individuals is NA, and should be ignored
          #when calculaing weighted mafs.
          dm[i,] <- 0 #likewise, reset the number of inds moving to zero to avoid messing up weights later
          #browser()
        }
        else{
          inds <- rbinom(nt, 2, maf[i]) #do a draws for each individual from the input maf
          pos <- 0 #set starting position in the inds vector to zero
          for(j in 1:(length(nds))){ #partition individuals form the inds vector into each location
            if(nds[j] == 0){bm[i,j] <- NA} #if no individuals disperesed here, set the allele count as NA and move on
            else{
              bm[i,j] <- sum(inds[(pos+1):(pos+nds[j])]) #get the allele count for each spot
              pos <- pos + nds[j]
            }
          }
        }
      }
      mafs <- colSums(bm, na.rm = TRUE)/(2*colSums(dm, na.rm = TRUE))
      #print(mafs)
      datm[,(h+2)] <- mafs
    }
    datm[,1]<-i.mat%*%datm[,1] #spread pop
    datm[,1]<-BH(datm[,1], r, k) #grow pop
  }
  colnames(datm) <- c("N", "xs", paste0("loci", 1:l)) #change column names for datm
  return(list(dat = as.data.frame(datm), imafs = imaf)) #return datm and imaf
}

#Matrix generating function for branch on branch. The kernal models two step movement from starting
#point (x) to the fork (y) and then to end point(z), with the average movement during each step
#the overall m multiplied by the relative distance traveled during both steps (always sums to m).
#takes a specific set of xs values, the stream positions which are on branches.
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
  return(out)
}
#fagen, william dendritic networks

#Same as A3, but takes a vector of y values (for multiple junctions). Still needs to be halved, quartered
#ect. for branch weights as necissary (should be preformed on resulting matrix).
A3 <- function(a,m,L,n,y,xs){
  K3 <- function(x,z){
    y <- c(x,y,z)
    part <- numeric(length(y) - 1)
    #print(y)
    for(i in 1:(length(y) - 1)){
      #print(i)
      trav <- abs(y[i] - y[i+1])
      s <- trav/sum(abs(y[-1] - y[-length(y)]))
      #cat("y[i], y [i+1]:", y[i], y[i+1], "\n")
      part[i] <- -abs(y[i+1] - (m*s) - y[i])
      #print(part)
    }
    out <- exp(sum(part)/a)/(2*a)
    return(out)
  }
  deltax <- 2*L/n #set range of x
  out <- outer(xs, xs, K3)*deltax #create matrix
  return(out)
}

#create.matrix
#Creates a disperal matrix based on the K3 kernal for a network of branches. Result is compatible with stream.drift when in.mat = TRUE.
#Inputs:
# brchs: A data frame with information on the branches in the following syntax:
#        <branch.id><branch.length><start.vertex><start.weight><end.vertex><end.weight>
#           branch.id: The name of the branch.
#           branch.length: The length of the branch.
#           start.vertex: ID for the upstream fork/confluence. Must be unique for each vertex. Terminal vertices must still be uniquely named.
#           start.weight: A numeric entry giving the weight of the branch at the upstream vertex relative to all other branches entering that vertex. Do not need to sum to one.
#           end.vertex: ID for the downstream fork/confluence. Must be unique for each vertex. Terminal vertices must still be uniquely named.
#           end.weight: A numeric entry giving the weight of the branch at the downstream vertex relative to all other branches entering that vertex. Do not need to sum to one.
#
# a: 2a^2 is the variance in displacement for dispersal
# m: displacement distance
# dx: delta x, the interval/step size between locations where dispersal is calibrated. 
#     The distances of all branches MUST be divisible by dx.
# give.xs: Logical. If TRUE, returns a two part list: 1) the dispersal matrix, 2) a data frame
#          where the first column is the position along each branch and the second is the branch ID.
#          These positions are the positions refered to in the dispersal matrix, in the order given.
# normalize: Logical. If TRUE, normalizes each column to sum to 1 (no individuals leave the system or enter the
#            system except by growth).

#Outputs:
# Dispersal matrix and (if give.xs = TRUE) a dataframe detailing the positions the dispersal matrix
# refers to.
#
#Notes:
# Undercase s and e are restricted vertex names, and will cause errors if used in input.
Edge.to.Matrix <- function(brchs, a, m, dx, give.xs = TRUE, normalize = TRUE){
  if(any(brchs[,2] %% dx != 0)){
    stop("All branch lengths MUST be divisible by the given dx value. If this is not the case,
        interpolation will not be consistant across branches. 0 length branches are acceptable, and
        will be counted as two vertices with weights.")
  }
  
  require(igraph)
  #First, designate functions to make branch on branch matrices.
  #A) No nodes
  #A3) Many nodes
  
  #simple matrix function for NO nodes
  A <- function(a, m, dx, xs){
    K <- function(y, x){#make kernal function
      out <- exp(-abs(y - m - x)/a)/(2*a)
      return(out)
    }
    out <- outer(xs, xs, K)*dx #create matrix
    return(out)
  }
  
  #Same as A, but takes a vector of y values (for multiple nodes). Still needs to be halved, quartered
  #ect. for branch weights as necissary (should be preformed on resulting matrix).
  A3 <- function(a,m,dx,y,xs.s,xs.e){
    K3 <- function(x,z){
      y <- c(x,y,z)
      part <- numeric(length(y) - 1)
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
      return(out)
    }
    out <- outer(xs.s, xs.e, Vectorize(K3))*dx #create matrix
    return(out)
  }
  
  #construct igraph
  cat("Constructing network...\n")
  bnet <- graph_from_edgelist(as.matrix(brchs[,c(3,5)]), directed = FALSE) #construct igraph network
  # V(bnet)$name <- 1:vcount(bnet) #name vertices
  if (give.xs == TRUE){ #if the xs data is requested...
    xs.output <- data.frame(xs = numeric(1), branch = character(1)) #initialize
    xs.output$branch <- as.character(xs.output$branch)
  }
  
  
  ##Next, loop through the branch input file and create matrix (large) row. After the first row, will append new row
  #loop through each row, do every combination (branch j influence on branch i)
  cat("Creating output matrix...\n")
  for(i in 1:nrow(brchs)){ #for each (ending) branch...
    if(brchs[i,2] == 0){ #skip to the next starting branch if no actual locations on this branch (too small)
      next
    }
    xs.e <- seq(dx/2, brchs[i,2] - dx/2, by = dx) #ending xs, most upstream at 0
    if(give.xs == TRUE){ #if the xs data is requested...
      xs.output <- rbind(xs.output, data.frame(xs = xs.e, branch = brchs[i,1])) 
    }
    mc <- matrix(NA, 1, length(xs.e)) #initialize output matrix
    for(j in 1:nrow(brchs)){ #loop through and compare to each (starting) branch...
      xs.e <- seq(dx/2, brchs[i,2] - dx/2, by = dx) #reset xs.e
      if(i == j){ #for self comparison
        mp <- A(a, m, dx, xs.e) #no nodes or weights for branch on itself
        mc <- rbind(mc, mp)
        next
      }
      if(brchs[j,2] == 0){ #skip to next ending branch if no actual locations on this branch (too small)
        next
      }
      else{ #for all others...
        #need to find all nodes passed through
        be.ns <- as.character(brchs[i,3]) #branch one start node
        be.ne <- as.character(brchs[i,5]) #branch one end node
        bs.ns <- as.character(brchs[j,3]) #branch two start node
        bs.ne <- as.character(brchs[j,5]) #branch two end node
        #print(bnet)
        #print(brchs[i,3])
        #cat(bs.ns, bs.ne, be.ns, be.ne, "\n")
        paths <- all_simple_paths(bnet + vertex(c("s", "e")) + 
          edges(c(bs.ns, "s"), c(bs.ne, "s"), c(be.ns, "e"), c(be.ne, "e")) -
          edges(get.edge.ids(bnet, c(bs.ns, bs.ne)), get.edge.ids(bnet, c(be.ns, be.ne))),
          "s", "e") #add temporary vertex midway through each edge, remove old edges to
                    #avoid loops, then get the simple paths from one new node (the starting
                    #point) to the other (the ending point).
        
        #next, need to construct the matrix for each path.
        mp <- list()
        for(k in 1:length(paths)){
          #for each path..
          #get the vertices
          verts <- paths[[k]][-c(1, length(paths[[k]]))] #get the true vertices, minus dummies
          
          #set initial xs
          xs.s <- seq(dx/2, brchs[j,2] - dx/2, by = dx) #caluclate xs for this branch (starting)
          
          #need to use A3 to generate matrix, which means getting "positions" for vertices
          #these are positions relative to starting and ending xs values. These will be - or +,
          #at the end of the branches.
          
          #get branch IDs between vertices
          t.edges <- get.edge.ids(bnet,c(rbind(V(bnet)$name[verts][-length(verts)], V(bnet)$name[verts][-1])), directed = FALSE)
          
          #now have edges passed through, and therefore lengths and direction
          #need to set y positions. Get weights at the same time
          y <- ifelse(brchs[j,3] == V(bnet)$name[verts[1]], 0, brchs[j,2]) #first y is zero if first element of verts is a "start" node, the length of the branch if it is an "end" node.
          w <- numeric(length(verts)) #initialize weights vector
          lvert <- ifelse(bs.ns == V(bnet)$name[verts[1]], bs.ne, bs.ns) #the initial "last vertex" is the vertex of the starting branch which isn't the first of the path vertices
          if(length(verts) > 1){
            for (h in 1:length(t.edges)){
              tvert <- names(verts)[h] #set the name of the the current vertex, save typing
              y[h+1] <- ifelse(brchs[t.edges[h],3] == tvert, y[h] + brchs[t.edges[h],2], y[h] - brchs[t.edges[h],2]) #if this vertex is a "start node", set y as the previous y plus the distance of the river, otherwise set it as y minus the distance of the river
              wi <- ifelse(brchs[t.edges[h],3] == tvert, brchs[t.edges[h],4], brchs[t.edges[h],6]) #if this vertex is and end node, use end weight, otherwise use start weight
              forks <- brchs[which((brchs[,3] == tvert | brchs[,5] == tvert) & 
                                  (brchs[,3] != lvert & brchs[,5] != lvert)),] #get the branches which don't involve the previous vertex but do involve the current one
              #print(forks)
              w[h] <- wi/sum(ifelse(forks[,3] == tvert, forks[,4], forks[,6])) #relative weight of the branch used in this path is its weight divided by the total weights. ifelse selects the correct weight (start or end) for the node.
              lvert <- tvert #set tvert to lvert for next loop
            }
          }
          else{
            tvert <- V(bnet)$name[verts[1]]
          }
          tvert <- V(bnet)$name[verts[length(verts)]]
          fvert <- ifelse(be.ns == V(bnet)$name[verts[length(verts)]], be.ne, be.ns) #set final destination vertex
          wi <- ifelse(brchs[j,3] == fvert, brchs[j,4], brchs[j,6]) #get the appropriate final weight
          forks <- brchs[which((brchs[,3] == tvert | brchs[,5] == tvert) & 
                                 (brchs[,3] != lvert & brchs[,5] != lvert)),] #get possible forks for current vertex, not involving previous vertex
          w[length(w)] <- wi/sum(ifelse(forks[,3] == tvert, forks[,4], forks[,6])) #get the weight of the final choice
          
          #change the ending xs 
          if(brchs[i,3] == tvert){ #if the final node in verts is a start node...
            xs.e <- xs.e + y[length(y)] #this is approaching from upstream, add the last value of y to xs.e
          }
          else{ #otherwise
            xs.e <- y[length(y)] - xs.e #this is approaching from downstream, take the last value of y and subtract each xs.e from it
            xs.e <- sort(xs.e)
          }
          
          #run A3 to get the matrix part
          mp[[k]] <- A3(a, m, dx, y, xs.s, xs.e)
          mp[[k]] <- mp[[k]]*prod(w) #weight the matrix by the total weights of each vertex fork
          #if(max(mp[[k]]) >= 1){browser()} #get this when a is greater than .5 and dispersal distance = m (including direction) 
        }
        #sum the different path matrixes to get final matrix
        mp <- Reduce("+", mp)
      }
      #need to position mp
      if(is.list(mp) == TRUE){
        mp <- unlist(mp)
      }
      mc <- rbind(mc, mp)
    }
    #need to position large row
    mc <- mc[-1,] #remove filler row column at start
    if(exists("full.matrix") == FALSE){
      full.matrix <- matrix(NA, nrow(mc), 1) #if it doesn't exist yet, intialize the large matrix, using the colnumbers of the first "row"
    }
    full.matrix <- cbind(full.matrix, mc) #bind the large matrix row to the full matrix
  }
  full.matrix <- full.matrix[,-1] #remove the filler column
  
  # normalize
  if(normalize){
    full.matrix <- full.matrix/rep(colSums(full.matrix), each = nrow(full.matrix))
  }
  
  if (give.xs == TRUE){
    xs.output <- xs.output[-1,] #cut off initialization row
    output <- list(full.matrix, xs.output) #combine the matrix and xs.output into a list
    return(output) #output both matrix and xs.output
  }
  else{
    return(full.matrix) #output the final matrix
  }
}



#function to create an edgefile for input into create.matrix from a SpatialLinesDataFrame.
#inputs: soi: SpatialLinesDataFrame containing streams.
#        ei: ElevationIndex, a SpatialLinesDataframe containing elevational contours for determining elevations of the vertices.
#        ei.c: Name or number of the column containing the value to index against in ei (elevation, ext). Higher values in this column will be designated as upstream.
#        plot.check: If TRUE, plots the soi streams and the vertices to check that they are correct.
GIS.to.Edge <- function(soi, ei, ei.c, plot.check = TRUE){
  require(sp)
  require(raster)
  require(rgdal)
  require(dplyr)
  require(ggplot2)
  
  
  cat("Reminder: after vertices are found, this function may require user inputs.\nDisaggregating...\n")
  #first need to fully disaggregate soi
  d.soi <- disaggregate(soi)
  while(!identical(soi, d.soi)){ #while not fully disaggregated...
    soi <- d.soi #reset soi to d.soi
    d.soi <- disaggregate(soi) #disaggretate again.
  }
  soi <- d.soi #reset once more
  remove(d.soi) #remove the extra copy
  
  
  
  #get coordinates for all of the possible vertices, both internal and external.
  #get the internal vertices:
  cat("Finding all vertices..\n")
  soi.c <- coordinates(soi) #get the coordinates of the lines
  wmat <- numeric(2)
  for(i in 1:length(soi.c)){ #rbind all the line coordinates together
    wmat <- rbind(wmat, do.call("rbind",soi.c[[i]]))
  }
  wmat <- wmat[-1,] #remove the intialization column
  colnames(wmat) <- c("x", "y")
  
  zd <- wmat[duplicated(wmat),] #get any points that appear more than once, internal vertices
  
  #get the start and end vertices as well:
  #from the sn function
  startEndPoints <- function(x) {
    firstLast <- function(L) {
      cc <- coordinates(L)[[1]]
      rbind(cc[1, ], cc[nrow(cc), ])
    }
    do.call(rbind, lapply(x@lines, firstLast))
  }
  
  ends <- startEndPoints(soi) #get the start and end points
  zd_es <- unique(rbind(zd, ends)) #rbind start and end points to internal vertices and get only unique vertices.
  rownames(zd_es) <- 1:nrow(zd_es)
  
  
  #figure out elevations for the vertices
  cat("Figuring out elevation for", nrow(zd_es), "total vertices...\n")
  temp <- SpatialPoints(zd_es, proj4string = CRS(proj4string(soi))) #make a sp df
  v.els <- numeric(nrow(data.frame(temp))) #initialize output vector
  branch_map_data <- ggplot2::fortify(soi)
  branch_map_data$vertex <- NA
  
  for (k in 1:length(v.els)) { #loop through the vertices. Could do this by grabbing two closest, interpolating if I wanted to be more accurate.
    cat("Vertex:", k, "\n")
    ds <- gDistance(temp[k,], ei, byid = TRUE) #get all of the distances to elements of the elevation SLDF
    ds <- cbind(1:length(ds), ds) #add the indices back
    colnames(ds) <- c("index", "dist") #add colnames for sorting
    ds <- as.data.frame(ds) 
    ds <- arrange(ds, dist) 
    ds <- ds[!duplicated(ds[,2]),] #take only unique distances, since repeats are almost certainly data errors. Could make an option.
    w1 <- 1-(ds[1,2]/(ds[1,2]+ds[2,2])) #get the weight based on proximity to the two closest points. Could do this for ALL points or trianglulate rather than the first two...
    w2 <- 1-w1
    wd <- w1*as.numeric(as.data.frame(ei[ds[1,1], ei.c])) + 
      w2*as.numeric(as.data.frame(ei[ds[2,1], ei.c])) #get the weighted distance
    v.els[k] <- wd #print the distance.
    
    
    # now, figure out which stream portions map to this vertex for later use
    matches <- which(apply(branch_map_data[,1:2], 1, function(x) all(x == zd_es[k,])))
    branch_map_data$vertex[matches] <- k
  }
  
  zd_es <- cbind(zd_es, v.els)
  colnames(zd_es)[3] <- "elev"
  
  #prepare zd_es data frame
  zd_es_df <- as.data.frame(zd_es)
  zd_es_df$ID <- 1:nrow(zd_es_df) #add ID
  zd_es_df$el.div <- zd_es_df$elev - mean(zd_es_df$elev) #get the elevation deviation
  off.scale.x <- max(zd_es[,1]) - min(zd_es[,1]) #set the x scale
  off.scale.y <- max(zd_es[,2]) - min(zd_es[,2]) #set the y scale
  zd_es_df$e.rank <- floor((rank(x = desc(zd_es_df$elev)))) #get the rank, from high to low elevation
  if(any(duplicated(zd_es_df$e.rank))){plot.check <- TRUE} #if there are any matches, MUST check the plot and give corrections
  
  while(plot.check){ #while plot.check is true, plot the intermediate nodes and lines

    #prepare a label
    zd_es_df$lab <- paste0(zd_es_df$ID, ",", zd_es_df$e.rank) #set the label to show, ID then rank
    
    #plot
    c.plot <- ggplot() + geom_path(data = fortify(soi), aes(x = long, y = lat, group = group), col = "grey") + 
      geom_point(data = zd_es_df, aes(x = x, y = y, color = el.div), size = 2) +
      scale_color_gradient2(high = "red", mid = "blue", low = "purple") +
      geom_text(data = zd_es_df, aes(x = x + (off.scale.x/90), y = y + (off.scale.y/90), label = lab), color = "darkred")
    
    print(c.plot)
    View(zd_es_df)
    cat("Vertex check plot and data.frame opened.\n")
    
    
    #accept inputs, rearrange point ranks manually.
    #check to see if adjustments are needed:
    if(any(duplicated(zd_es_df$e.rank))){
      duplist <- zd_es_df[zd_es_df$e.rank %in% unique(zd_es_df$e.rank[duplicated(zd_es_df$e.rank)]),]$lab
      cat("Points with equal elevation ranks detected:", duplist, "Please Fix:\n")
      resp  <- "y"
    }
    else{
      cat("Check Plot: Vertex labels are vertex ID followed by vertex elevation rank.\nAre adjustments needed? (y or n)")
      resp <- readLines(n = 1)
    }
    if(resp == "y"){
      cat("Specify rearrangements in format <ID,ID><space><ID,ID,ID><space><ID,ID>... on one line, where the upstream ID is first and the downstream ID is second.")
      resp <- readLines(n = 1) #get the input
      resp <- unlist(strsplit(resp, " ")) #split the input
      for(l in 1:length(resp)){
        dup.sets <- FALSE #reset dup.sets for later logical
        t.or <- unlist(strsplit(resp[l], ",")) #split this reorder
        t.points <- zd_es_df[zd_es_df$ID %in% t.or,] #get points to rearrange
        other.points <-  zd_es_df[!(zd_es_df$ID %in% t.or),] #get all of the points not in this set

        if(any(duplicated(t.points$e.rank))){#if there are duplicate points, grab the duplicate sets, fix them in the order desired.
            dups <- t.points[t.points$e.rank %in% unique(t.points$e.rank[duplicated(t.points$e.rank)]),] #get any duplicated points only
            dup.sets <- unique(dups$e.rank) #elevation rank options that have duplicates to fix
            fixed.dups <- list()
            for(k in 1:length(dup.sets)){ #for each dup set
              t.set <- dups[dups$e.rank == dup.sets[k]]
              given.order <- t.or[which(t.or %in% t.set$ID)]
              r.set <- numeric(ncol(t.set)) #initialize
              for(j in 1:length(given.order)){ #for each element in this...
                r.set <- rbind(r.set, t.set[t.set$ID == given.order[j],]) #add each reordered element sequentially
              }
              r.set <- r.set[-1,] #remove filler row
              r.set$new.ranks <- seq(1, nrow(r.set), 1) #re-designate ranks for within family
              #next, save fixed dups
              fixed.dups[[paste0("dup", dup.sets[k])]] <- r.set #add the set to the list of fixed dfs, named by the e.rank which was duplicated.
            }
        }
        
        #now, rearrange UNIQUE (not the fixed dups) points according to order given, then add in the duplicate fixed points.
        #don't add a new.rank to it yet
        u.t.points <- t.points[!duplicated(t.points$e.rank),] #unique points.
        u.t.or <- t.or[t.or %in% u.t.points$ID] #order of the unique points to use
        u.t.or <- data.frame(cbind(u.t.or, seq(1:length(u.t.or)))) #get the order explicitly tied to the ID
        colnames(u.t.or) <- c("ID", "sortby") 
        u.t.points <- merge.default(u.t.points, u.t.or, by = "ID") #merge the order info with the points to sort
        u.t.points <- arrange(u.t.points, sortby) #arrange the points by the order to sort
        u.t.points <- u.t.points[,-ncol(u.t.points)]
        u.t.points$new.ranks <- sort(u.t.points$e.rank) #resort to e.ranks into the new order, save as new.ranks
        
        #now resorted. Now just need to paste everything together
        
        fixed.order <- numeric(ncol(zd_es_df)) #initialize
        unique.points <- nrow(u.t.points)
        for(k in 1:unique.points){
          #add anything less than this point 
          fixed.order <- rbind(fixed.order, other.points[other.points$e.rank < min(u.t.points$new.ranks),]) #add anything less than this point.
          other.points <- other.points[other.points$e.rank > min(u.t.points$new.ranks),] #keep only the rest
          
          #check to see if the point to add is part of a duplicate family, if so, add the duplicate family and increase everything else's rank
          if(dup.sets){
            if (any(u.t.points$e.rank[k] %in% dup.sets)){ #is this rank in a duplicate family?
              t.part <- fixed.dups[[paste0("dup",u.t.points$e.rank[k])]] #rbind this fixed dup family
              t.part$e.rank <- t.part$new.ranks #set e.rank to new ranks
              t.part <- t.part[,-ncol(t.part)] #remove the new ranks column
              t.part$e.rank <- t.part$e.rank + max(fixed.order$e.rank) #change the rank to be the family rank plus the rank already down.
              add.rank <- nrow(t.part) - 1 #get the number added
              fixed.order <- rbind(fixed.order, t.part) #add this part.
              other.points$e.rank <- other.points$e.rank + add.rank #add the additional ranks to the other points 
              u.t.points$new.ranks <- u.t.points$new.ranks + add.rank #add the additional ranks to the remaining u.t.points
            }
          }
          else{
            t.part <- u.t.points[1,] #get this part, the first in each case since we are chopping off
            t.part$e.rank <- t.part$new.ranks #rename it 
            t.part <- t.part[,-ncol(t.part)] #remove old new.ranks column
            fixed.order <- rbind(fixed.order, t.part) #add this part
          }
          #chop off this point from u.t.points
          u.t.points <- u.t.points[-1,] #remove the row we just delt with.
        }
        #add the remaining other points
        fixed.order <- rbind(fixed.order, other.points)
        fixed.order <- fixed.order[-1,] #remove the initialization
        zd_es_df <- fixed.order #reset zd_es_df with the new ordering and go to the next reorder.
      }
    }
    else{plot.check <- FALSE} #if no re-arrangements called for, set plot.check to false and move on
  }
  
  
  
  #now, need to split the streams at these points, get lengths, add to output
  cat("Splitting streams at vertices and outputing edge info...\n")
  out <- data.frame(branch.id = NA, branch.length = NA, start.vertex = NA,
                    start.weight = NA, end.vertex = NA, end.weight = NA)
  s.c <- 1 #initialize stream count
  for(i in 1:length(soi.c)){ #for each stream portion
    cat("Working on stream", i, "of", length(soi.c), "\n")
    ms <- numeric(2) #initialize matching streams info
    
    #determine if any coordinates in this stream section match any of the starting/ending/internal vertices. Save info if so.
    for(j in 1:nrow(zd_es)){ #for each vertex
      xmatch <- which(soi.c[[i]][[1]][,1] == zd_es[j,1]) #get any x.coord matches for the vertex
      if(length(xmatch) != 0){ #if there were any matches... 
        for (k in 1:length(xmatch)){ #for each match...
          if(soi.c[[i]][[1]][xmatch[k],2] == zd_es[j,2]){ #does the y.coord also match?
            ms <- rbind(ms, cbind(xmatch[k], j)) #if so, save the index for the stream and the list of vertices
          }
        }
      }
    }
    ms <- ms[-1,] #remove the initialization row
    if(nrow(ms) <= 1){stop} #error, since there should always be at least two matches.
    
    #arrange the vertices in the stream segement by the index in the stream.
    colnames(ms) <- c("s.index", "v.index")
    ms <- as.data.frame(ms)
    ms <- arrange(ms, s.index) #arrange by stream index. this will sort them and allow for easy stream splitting
    
    
    #split the stream up at each vertex, save info for export.
    for(j in 1:(nrow(ms)-1)){#for each match
      segment <- soi.c[[i]][[1]][ms[j,1]:ms[(j + 1),1],] #this gets the portion of the stream between this vertex and the next
      out[s.c, "branch.id"] <- s.c
      out[s.c, "branch.length"] <- LineLength(segment)
      out[s.c, "start.vertex"] <- ifelse(zd_es_df[zd_es_df$ID == ms[j,2],6] < zd_es_df[zd_es_df$ID == ms[j+1,2],6], ms[j,2], ms[j+1,2]) #print out higher elevation vertex
      out[s.c, "end.vertex"] <- ifelse(zd_es_df[zd_es_df$ID == ms[j,2],6] > zd_es_df[zd_es_df$ID == ms[j+1,2],6], ms[j,2], ms[j+1,2]) #print out lower elevation vertex
      s.c <- s.c + 1 #increase the stream count
    }
  }
  
  # clean up the branch map data, work point
  cbmd <- na.omit(branch_map_data)
  branch_map_data$branch <- 0
  for(i in 1:(nrow(cbmd) - 1)){
    if(as.numeric(rownames(cbmd)[i]) + 1 == as.numeric(rownames(cbmd)[i+1])){next} # if we are moving on to the next segment (aka, no stream segments between this vertex and the next), skip
    match.branch <- c(cbmd[i,]$vertex, cbmd[i + 1,]$vertex)
    match.branch <- which(out$start.vertex %in% match.branch & out$end.vertex %in% match.branch)
    if(length(match.branch) == 0){browser()}
    branch.rows <- as.numeric(row.names(cbmd)[i:(i+1)])
    branch_map_data$branch[branch.rows[1]:branch.rows[2]] <- match.branch
  }
  
  # convert the matrix IDs to alpha
  out$start.vertex <- paste0("vertex_", out$start.vertex)
  out$end.vertex <- paste0("vertex_", out$end.vertex)
  branch_map_data$vertex[-which(is.na(branch_map_data$vertex))] <- paste0("vertex_", branch_map_data$vertex[-which(is.na(branch_map_data$vertex))])
  
  return(list(edges = out, plot = c.plot, map_data = branch_map_data))
}



