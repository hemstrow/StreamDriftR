library(sp)
library(igraph)

setClass("igraph")
setClass("SpatialLinesNetwork", representation(sl = "SpatialLines", g = "igraph", 
                                               nb = "list"), validity = function(object) {
                                                 stopifnot(length(object@sl) == length(E(object@g)))
                                                 stopifnot(length(object@nb) == length(V(object@g)))
                                               })


SpatialLinesNetwork = function(sl) {
  stopifnot(is(sl, "SpatialLines"))
  if (!is(sl, "SpatialLinesDataFrame")) 
    sl = new("SpatialLinesDataFrame", sl, data = data.frame(id = 1:length(sl)))
  if (!all(sapply(sl@lines, length) == 1)) 
    stop("SpatialLines is not simple: each Lines element should have only a single Line")
  #need to get the different nodes
  sl <- disaggregate(sl)
  soi.c <- coordinates(sl) #get the coordinates of the lines
  wmat <- c(NA,NA)
  for(i in 1:length(soi.c)){ #rbind all the line coordinates together
    wmat <- rbind(wmat, do.call("rbind",soi.c[[i]]))
  }
  wmat <- wmat[-1,]
  zd <- wmat[duplicated(wmat),]
  #this works, still need to get end points
  
  #from the sn function
  startEndPoints <- function(x) {
    firstLast <- function(L) {
      cc <- coordinates(L)[[1]]
      rbind(cc[1, ], cc[nrow(cc), ])
    }
    do.call(rbind, lapply(x@lines, firstLast))
  }
  
  
  s <- startEndPoints(sl) #get the start and end points
  zd <- unique(rbind(zd, s)) #rbind this to the nodes and get only the unique start and end points
  
  pts = 1:nrow(s)
  
  # the following can't be done vector-wise, there is a progressive effect:
  if (nrow(zd) > 0) {
    for (i in 1:nrow(zd)) pts[zd[i, 2]] = pts[zd[i, 1]]
  }
  stopifnot(identical(s, s[pts, ]))
  
  # map to 1:length(unique(pts))
  pts0 = match(pts, unique(pts))
  node = rep(1:length(sl), each = 2)
  nb = lapply(1:length(unique(pts)), function(x) node[which(pts0 == x)])
  g = graph(pts0, directed = FALSE)  # edges
  nodes = s[unique(pts), ]
  g$x = nodes[, 1]  # x-coordinate vertex
  g$y = nodes[, 2]  # y-coordinate vertex
  g$n = as.vector(table(pts0))  # nr of edges
  # line lengths:
  sl$length = sapply(sl@lines, function(x) LineLength(x@Lines[[1]]))
  E(g)$weight = sl$length
  # create list with vertices, starting/stopping for each edge?  add for
  # each SpatialLines, the start and stop vertex
  pts2 = matrix(pts0, ncol = 2, byrow = TRUE)
  sl$start = pts2[, 1]
  sl$end = pts2[, 2]
  new("SpatialLinesNetwork", sl = sl, g = g, nb = nb)
}
