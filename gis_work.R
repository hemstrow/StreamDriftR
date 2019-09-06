library(sp)
library(rgdal)
library(raster)
library(RColorBrewer)
library(rgeos)
library(maptools)
source("R/sn2.R"); source("R/stream_drift.R")

# #stream and watershed infiles
# rivers <- shapefile("shapefiles/rivers/rivers.shp")
# ws <- shapefile("shapefiles/watersheds/Oregon_Watershed_Councils_2014.shp")
# dcws <- ws[ws$altName == "Upper Deschutes WC" | ws$altName == "Crooked River WC",]
# #plot(dcws)
# dcws <- aggregate(dcws)
# newCRS <- CRS("+proj=lcc +lat_1=43 +lat_2=45.5 +lat_0=41.75 +lon_0=-120.5 +x_0=399999.9999984001 +y_0=0 +datum=NAD83 +units=ft +no_defs +ellps=GRS80 +towgs84=0,0,0")
# 
# #elevation infiles, crop them, fix their row ids, and rbind them.
# e_c <- shapefile("shapefiles/contours/Contour_100ce/Contour_100ce.shp")
# e_s <- shapefile("shapefiles/contours/Contour_100s/Contour_100s.shp")
# e_w <- shapefile("shapefiles/contours/Contour_100w/Contour_100w.shp")
# e_c_c <- crop(e_c, dcws)
# e_s_c <- crop(e_s, dcws)
# e_w_c <- crop(e_w, dcws)
# e_c_c <- spChFIDs(e_c_c, as.character(1:nrow(e_c_c)))
# at.row <- nrow(e_c_c)
# e_s_c <- spChFIDs(e_s_c, as.character((1+at.row):(at.row + nrow(e_s_c))))
# at.row <- at.row + nrow(e_s_c)
# e_w_c <- spChFIDs(e_w_c, as.character((1+at.row):(at.row + nrow(e_w_c))))
# elev <- spRbind(e_c_c, e_s_c)
# elev <- spRbind(elev, e_w_c)
# elev <- spTransform(elev, newCRS)
# remove(e_c, e_s, e_w, e_c_c, e_s_c, e_w_c)
# 
# 
# #crop rivers
# rivers <- spTransform(rivers, newCRS)
# dr <- crop(rivers, dcws)
# remove(rivers)
# 
# # plot(dr, col = "red")
# # plot(dcws, add = TRUE)
# # plot(soi, add = TRUE)
# # plot(elev[as.numeric(elev$CONTOUR) == 2000 | as.numeric(elev$CONTOUR) == 2100,], add = TRUE, col = "grey")
# 
# # plot(dr)
# # plot(dcws, add = TRUE)
# 
# soi <- dr[dr$NAME == "Pilot Butte Canal" |
#             dr$NAME == "Deschutes River" |
#             dr$NAME == "Crooked River"|
#             dr$STREAMS_ID == "18191",]
# soi <- disaggregate(soi)
# 
# # plot(soi, col = "red")
# # plot(dcws, add = TRUE)
# # plot(elev[as.numeric(elev$CONTOUR) %% 500 == 0,], col = brewer.pal(9, "Greys"), add = TRUE)
# 
# map.dat <- list(soi = soi, ei = elev, dr = dr)
# saveRDS(map.dat, "shapefiles/map_dat.RDS")



#=======================start here on a second approach===================
map.dat <- readRDS("shapefiles/map_dat.RDS")

# get edges
deschutes <- GIS.to.Edge(map.dat$soi, ei = map.dat$ei, ei.c = "CONTOUR", TRUE)
d.plot <- deschutes$plot
d.map.dat <- deschutes$map_data
deschutes <- deschutes$edges

deschutes$branch.length <- deschutes$branch.length/5200
deschutes$start.weight <- 1
deschutes$end.weight <- 1
deschutes$branch.length <- round(deschutes$branch.length)

# make the matrix
deschutes.m <- Edge.to.Matrix(deschutes, a = 2, m = 5, dx = 1)

# initial pops, 100 at top of system
n0 <- numeric(nrow(deschutes.m[[2]]))
n0[which(deschutes.m[[2]]$branch == 7)][1] <- 100

# run
d.drift <- stream.drift(n0, l = 100, r = 2, k = 500, Tf = 100, in.mat = deschutes.m[[1]], in.xs = deschutes.m[[2]]$xs)

out <- d.drift
# add branch metadata back in
out$dat$branch <- deschutes.m[[2]]$branch

# melt and plot
mdat <- reshape2::melt(out$dat, id.vars = c("N", "xs", "branch"))
colnames(mdat)[c(4,5)] <- c("loci", "maf")
mdat$delta <- abs(mdat$maf - out$imafs[as.numeric(substr(mdat$loci, 5, 5))]) 
ggplot(mdat, aes(x = xs, y = delta, color = loci)) + geom_line() + theme_bw() + facet_wrap(~branch, scales = "free_x")

ave.diff <- tapply(mdat$delta, mdat[,2:3], mean, na.rm = T)
ave.diff <- reshape2::melt(ave.diff)
ave.diff <- na.omit(ave.diff)
ggplot(ave.diff, aes(x = xs, y = value, color = as.factor(branch))) + geom_line() + 
  scale_color_viridis_d() +
  directlabels::geom_dl(aes(label = as.factor(branch)), 
                        method = list(directlabels::dl.combine("first.points", "last.points"), cex = 0.8))



brchs <- unique(ave.diff$branch)
d.map.dat$delta
for(i in 1:length(brchs)){
  t.dat <- ave.diff[ave.diff$branch == brchs[i],]
  t.m.d <- d.map.dat[d.map.dat$branch == brchs[i],]
  sv <- t.m.d[1,]$vertex
  
  
  points.per.xs <- nrow(t.m.d)/nrow(t.dat)
  vals <- rep(t.dat$value, each = ceiling(points.per.xs))
  if(length(vals) > nrow(t.m.d)){
    rm.points <- ceiling(seq(1, length(vals), length = length(vals) - nrow(t.m.d)))
    if(rm.points[length(rm.points)] > length(vals)){
      rm.points[length(rm.points)] <- length(vals)
    }
    vals <- vals[-rm.points]
  }
  
  
  if(deschutes[deschutes$branch.id == brchs[i], "start.vertex"] == sv){
    d.map.dat[d.map.dat$branch == brchs[i],]$delta <- vals
  }
  else if(deschutes[deschutes$branch.id == brchs[i], "end.vertex"] == sv){
    d.map.dat[d.map.dat$branch == brchs[i],]$delta <- rev(vals)
  }
  else{
    stop("No vertex match.\n")
  }
}


ggplot(d.map.dat, aes(x = long, y = lat, color = delta, group = group)) + geom_path() + theme_bw() + scale_color_viridis_c()

#barriers <- shapefile("ODFW_44_5_ofpbds_shp/ofpbds_pt.shp")
#barriers <- spTransform(barriers, newCRS)
#barriers <- crop(barriers, dcws)
#impass <- barriers[barriers$fpbFPasSta == "Blocked"|barriers$fpbFPasSta == "Partial",]
#plot(impass)
#plot(soi, add = TRUE)

#wbs <- shapefile("waterbodies/waterbodies.shp")
#wbs <- spTransform(wbs, newCRS)
#plot(soi)
#plot(impass, add = TRUE)
#plot(wbs, add = TRUE)

#sn <- SpatialLinesNetwork(soi, 20000)

#plot(sn@g)
#plot(sn@g$x, sn@g$y, col = sn@g$n, pch = 16, cex = 2, asp = 1)
#lines(soi)
#text(sn@g$x, sn@g$y, E(sn@g), pos = 4)
