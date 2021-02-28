library(raster)
library(geosphere)
soi<- map_dat$soi
data<- raster("~/Documents/github/Fish/shapefiles/rasters.gri")
#barriers<- rgdal::readOGR("shapefiles/barriers/ofpbds_gdb.gdb", "ofpbds_pt")

river<- readRDS("shapefiles/transformed_river_spdf.RDS")
crs(data)<- crs(river)
#sf::st_as_sf(river)
#plot(sf::st_as_sf(river))
#segmentizeriver<- sf::st_segmentize(sf::st_as_sf(river), .1)

#sfbarriers<- sf::st_transform(sf::st_as_sf(barriers), as.character(crsriver))
#barriers<- sf::as_Spatial(sfbarriers)
#b<- crop(barriers, river)
river@lines[[1]]@Lines[[1]]@coords
#river.name<- river$NAME == "Pilot Butte Canal"
length<- .1
#n.parts<- 200

###### use st_snap or st_intersection for finding which barriers are where

#function to create an edgefile for input into create.matrix from a SpatialLinesDataFrame.
#inputs: soi: SpatialLinesDataFrame containing streams.
#        plot.check: If TRUE, plots the soi streams and the vertices to check that they are correct.
#        data: raster of geographic area where streams are located
#        river: SpatialLinesDataFrame of river
#        river.name: name of river
#        length: numerical value by which to divide soi (a larger number gives you less stream segments, smaller gives you more)

GIS.to.Edge <- function(soi, length, data, river, river.name, plot.check = TRUE){
  
  #grid?
  #alt: track distances of segments and put those into adj matrix
  
  #first we find the nodes
  #split streams into evenly-spaced segments where the fish live
  #make segments: http://rstudio-pubs-static.s3.amazonaws.com/10685_1f7266d60db7432486517a111c76ac8b.html
  #MergeLast <- function(lst) { #instead of merge last, divide last segment into previous ones
      l <- length(lst)
      lst[[l - 1]] <- rbind(lst[[l - 1]], lst[[l]])
      lst <- lst[1:(l - 1)]
      return(lst)
    }
  CreateSegments <- function(coords, length = 0, n.parts = 0) {
      #install geosphere package- create vector of lengths in segmentspatial lines function
      stopifnot((length > 0 || n.parts > 0))#get rid of this
      # calculate total length line
      total_length <- 0
      for (i in 1:(nrow(coords) - 1)) {
        d <- sqrt((coords[i, 1] - coords[i + 1, 1])^2 + (coords[i, 2] - coords[i + 
                                                                                 1, 2])^2)
        total_length <- total_length + d
      }
      
      # calculate stationing of segments
      
      #######ask will: should we set a check to make sure length is small enough to have make sure to have enough segments to make the fractions of the remainders small enough but still have genetic variance
      ###### check length to be half the length of smallest stream segment, if no length selected, set the length to half of smallest river segment or smaller
      #### sanity checks
      #do a default length - divide geosphere::lengthLine by 2
      if (length > 0) { #insert a check to make sure length is >0 but <whole segment #get rid of the remainders by putting them into each segment
        stationing <- c(seq(from = 0, to = total_length, by = length), total_length)
        r<- stationing[length(stationing)] - stationing[length(stationing)-1]
        newstationing<- stationing[1:(length(stationing)-1)]
        f<- r/(length(newstationing)-1)
        newstationing<- seq(0, f*(length(newstationing)-1), by = f) + newstationing
        
        #remainder (r) = last entry of list - second to last entry of list
        r <- stationing[length(stationing)] - stationing[length(stationing)-1]
        #make a new list (l) and input every single entry except for last entry 
        newstationing <- stationing[1:(length(stationing)-1)]
        #fraction (f) = r divided by number of entries of l-1
        f<- r/(length(newstationing)-1)
        #add f to each entry of l
        newstationing<- seq(0, f*(length(newstationing) - 1), by = f) + newstationing
        #fraction (f) = r divided by number of entries of l
        #add f to each entry of l

      } else {#get rid of this- unnecessary
        stationing <- c(seq(from = 0, to = total_length, length.out = n.parts), #this has the river divided into a number of parts, which is problem bc our segments are diff lengths
                        total_length)
      }
      
      # calculate segments and store the in list
      newlines <- list()
      for (i in 1:(length(stationing) - 1)) {
        newlines[[i]] <- CreateSegment(coords, stationing[i], stationing[i + 
                                                                           1])
      }
      return(newlines)
    }
  CreateSegment <- function(coords, from, to) {
      distance <- 0
      coordsOut <- c()
      biggerThanFrom <- F
      for (i in 1:(nrow(coords) - 1)) {
        d <- sqrt((coords[i, 1] - coords[i + 1, 1])^2 + (coords[i, 2] - coords[i + 
                                                                                 1, 2])^2)  #get rid of this (already done)
        distance <- distance + d
        if (!biggerThanFrom && (distance > from)) {#probably don't want it to do this but not a big deal
          w <- 1 - (distance - from)/d
          x <- coords[i, 1] + w * (coords[i + 1, 1] - coords[i, 1])
          y <- coords[i, 2] + w * (coords[i + 1, 2] - coords[i, 2])
          coordsOut <- rbind(coordsOut, c(x, y))
          biggerThanFrom <- T
        }
        if (biggerThanFrom) {
          if (distance > to) { #probably don't want it to do this but not a big deal
            w <- 1 - (distance - to)/d
            x <- coords[i, 1] + w * (coords[i + 1, 1] - coords[i, 1])
            y <- coords[i, 2] + w * (coords[i + 1, 2] - coords[i, 2])
            coordsOut <- rbind(coordsOut, c(x, y))
            break
          }
          coordsOut <- rbind(coordsOut, c(coords[i + 1, 1], coords[i + 1, 
                                                                   2]))
        }
      }
      return(coordsOut)
  }
  ### write sanity checks in here as well as Geosphere::lengthLine function
  
  
  SegmentSpatialLines <- function(sl, length = 0, n.parts = 0, merge.last = FALSE) {
      #put that down with the default length part
      stopifnot((length > 0 || n.parts > 0))
    
    ######test sanity check
    msg <- character()
    
    #length must be greater than 0/value
    if(length <= 0){
      msg <- c(msg, "Length must be greater than 0/value")
    }
    
    #n.parts must be greater than 0/value
    if(n.parts <= 0){
      msg <- c(msg, "n.parts must be greater than 0/value")
    }
    
    if(length(msg) > 0){
      stop(paste0(msg, collapse = "\n"))
    }
    
    return("All values are valid")
    #### end sanity check
    
      id <- 0
      newlines <- list()
      sl <- as(sl, "SpatialLines")
      #find lengths of all lines
      #for loop?
      lineLengths<- vector("numeric",length(sl))
      for(i in 1:length(sl)){
        lineLengths[i]<- geosphere::lengthLine(sl@lines[[i]]@Lines[[1]]@coords)#divide by 1000
      }
      defaultLength<- min(lineLengths)/2
      if(length > defaultLength){
        warning("Length must be less than or equal to half of the length of smallest line (in km). Using default length")
        length<- defaultLength
      }
      #if  length is too long, use default
      #for example if segment length is 2.1 km, and length provided is 
      for (lines in sl@lines) {
        for (line in lines@Lines) {
          crds <- line@coords
          # create segments
          segments <- CreateSegments(coords = crds, length, n.parts)
          #get rid of this
          if (merge.last && length(segments) > 1) {
            # in case there is only one segment, merging would result into error
            segments <- MergeLast(segments)
          }
          # transform segments to lineslist for SpatialLines object
          for (segment in segments) {
            newlines <- c(newlines, Lines(list(Line(unlist(segment))), ID = as.character(id)))
            id <- id + 1
          }
        }
      }
      return(SpatialLines(newlines))
    }
    
  
  #segments
  spdf<- SegmentSpatialLines(river, length, merge.last = TRUE)
  plot(spdf, col = rep(c("green", "red"), length.out = length(spdf)), axes = T)
  #segment coordinates
  spdf.c <- sp::coordinates(spdf)
  #remove duplicated points from spdf- use "duplicated", loop through spdf.c, spdf.c[-which(duplicated(spdf.c)),]
  for(i in length(spdf.c)){
    spdf.c[[i]][[1]]<- spdf.c[[i]][[1]][-which(duplicated(spdf.c)),]
  }
  
  wmat <- numeric(2)
  for(i in 1:length(spdf.c)){ #rbind all the line coordinates together
    wmat <- rbind(wmat, do.call("rbind",spdf.c[[i]]))
  }
  wmat <- wmat[-1,] #remove the intialization column
  colnames(wmat) <- c("x", "y")
  
  zd <- wmat[duplicated(wmat),] #get any points that appear more than once, internal vertices
  startEndPoints <- function(x) {
    firstLast <- function(L) {
      cc <- sp::coordinates(L)[[1]]
      rbind(cc[1, ], cc[nrow(cc), ])
    }
    do.call(rbind, lapply(x@lines, firstLast))
  }
  
  ends <- startEndPoints(spdf) #get the start and end points
  zd_es <- unique(rbind(zd, ends)) #rbind start and end points to internal vertices and get only unique vertices.
  rownames(zd_es) <- 1:nrow(zd_es)
  
  #next we figure out elevations
  #get elev data on river from "data" raster
  #slowest part so far- takes 17680 ms
  #try to just get the elevation points that match with vertices
  #maybe make a data frame from the vertices
  
  #try to parallelize
  elev.data<- raster::extract(data, zd_es)
  

  zd_es <- cbind(zd_es, elev.data)
  colnames(zd_es)[3] <- "elev"
  
  

  
  #prepare zd_es data frame
  zd_es_df <- as.data.frame(zd_es)
  zd_es_df$ID <- 1:nrow(zd_es_df) #add ID
  zd_es_df$el.div <- zd_es_df$elev - mean(zd_es_df$elev) #get the elevation deviation
  off.scale.x <- max(zd_es[,1]) - min(zd_es[,1]) #set the x scale
  off.scale.y <- max(zd_es[,2]) - min(zd_es[,2]) #set the y scale
  zd_es_df$e.rank <- floor((rank(x = dplyr::desc(zd_es_df$elev)))) #get the rank, from high to low elevation
  if(any(duplicated(zd_es_df$e.rank))){plot.check <- TRUE} #if there are any matches, MUST check the plot and give corrections
  
  #MAKE SURE HOW TO TELL WHICH NODES ARE ADJACENT TO EACH OTHER (nodes are the vertices in zd_es)
  #pseudocode
  #could use adjacent function in raster to tell which cells are adjacent
  #adjacent(data, )
  #then figure out which coordinates from spdf.c are in which cell
  #make a list of which vertices are in whihc stream segment
  # % in %
  #if true, return zd_es row and segment number of spdf.c
  #for looop??????
  #list with 4 columns (1st 2 vertex, 3rd is elev, 4th is stream segment)
  
  #should I call these edges?
  adjacent_nodes<- vector("list", length(spdf.c))
  for(i in 1:length(spdf.c)){
    adjacent_nodes[[i]]<- which(zd_es_df[,1] %in% spdf.c[[i]][[1]] & zd_es_df[,2] %in% spdf.c[[i]][[1]])
  }
  
  #####REPLACE:
  ####figure out how to make sure that the elevations go from high to low
  #which is the highest elev point (should be at some end of the segment)
  #psuedocode
  #get the adjacent nodes from above
  #use zd_es_df to see the elevation point corresponding to the node
  #if(!  elev.data point for node >  elev.data point for next node >  elev.data point for next node), 
  #send warning to user and have them fix
  #else move on
  #check three points
  #use for loop
  
  #for (each segment)
  #return (vertex with highest elevation )
  
  #for i in 1:length(adjacent_nodes)
  #for j in 1:length(adjacent_nodes[[i]])
  #
  
  #find adjacent segments
  #find average elevations on each segment
  #find diff in elev btwn segs
  #put into adjacency matrix
  #connectivity = barriers (like waterfalls)
  
  avg_elev <- vector("list", length(adjacent_nodes))
  diff_elev <- vector("list", length(adjacent_nodes))
  for(i in 1:length(adjacent_nodes)){
    #find average elev of each segment
    avg_elev[[i]]<- (sum(zd_es_df[adjacent_nodes[[i]], 3]))/(length(adjacent_nodes[[i]]))
    #find diff in elev btwn segs
    diff_elev[[i]]<- avg_elev[[i+1]] - avg_elev[[i]]
  }
  
  #for (each segment)
  #if (elev)
  #confused about whether i'm supposed to check the elevations of vertices in each segment or between two segments.i feel like btwn two segments doesn't really make sense
  
  
  while(plot.check){ #while plot.check is true, plot the intermediate nodes and lines
    
    #prepare a label
    zd_es_df$lab <- paste0(zd_es_df$ID, ",", zd_es_df$e.rank) #set the label to show, ID then rank
    
    #plot
    c.plot <- ggplot2::ggplot() + 
      ggplot2::geom_path(data = ggplot2::fortify(soi), ggplot2::aes(x = long, y = lat, group = group), col = "grey") + 
      ggplot2::geom_point(data = zd_es_df, ggplot2::aes(x = x, y = y, color = el.div), size = 2) +
      ggplot2::scale_color_gradient2(high = "red", mid = "blue", low = "purple") +
      ggplot2::geom_text(data = zd_es_df, ggplot2::aes(x = x + (off.scale.x/90), y = y + (off.scale.y/90), label = lab), color = "darkred")
    
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
        u.t.points <- dplyr::arrange(u.t.points, sortby) #arrange the points by the order to sort
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
          u.t.points <- u.t.points[-1,] #remove the row we just dealt with.
        }
        #add the remaining other points
        fixed.order <- rbind(fixed.order, other.points)
        fixed.order <- fixed.order[-1,] #remove the initialization
        zd_es_df <- fixed.order #reset zd_es_df with the new ordering and go to the next reorder.
      }
    }
    else{plot.check <- FALSE} #if no re-arrangements called for, set plot.check to false and move on
  }
  
 
  #rewrite?
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
    ms <- dplyr::arrange(ms, s.index) #arrange by stream index. this will sort them and allow for easy stream splitting
    
    
    #split the stream up at each vertex, save info for export.
    for(j in 1:(nrow(ms)-1)){#for each match
      segment <- soi.c[[i]][[1]][ms[j,1]:ms[(j + 1),1],] #this gets the portion of the stream between this vertex and the next
      out[s.c, "branch.id"] <- s.c
      out[s.c, "branch.length"] <- sp::LineLength(segment)
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

