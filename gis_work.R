library(sp)
library(rgdal)
library(raster)
#library(RColorBrewer)
library(rgeos)
# library(maptools)
source("R/sn2.R"); source("R/stream_drift.R")

# #stream and watershed infiles
rivers <- shapefile("shapefiles/rivers/rivers.shp")
ws <- shapefile("shapefiles/watersheds/Oregon_Watershed_Councils_2014.shp")
dcws <- ws[ws$altName == "Upper Deschutes WC" | ws$altName == "Crooked River WC",]
plot(dcws)
dcws <- aggregate(dcws)
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
rivers <- spTransform(rivers, crs)
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

# get and clean edges
deschutes <- GIS.to.Edge(map.dat$soi, ei = map.dat$ei, ei.c = "CONTOUR", TRUE) # Needed re-arrangements: 5,14 6,12 10,13
deschutes$edges$branch.length <- deschutes$edges$branch.length/5200
deschutes$edges$start.weight <- 1
deschutes$edges$end.weight <- 1
deschutes$edges$branch.length <- round(deschutes$edges$branch.length)

# make the matrix
deschutes.m <- Edge.to.Matrix(deschutes$edges, a = 2, m = 5, dx = 1)

# initial pops, 100 at top of system
n0 <- numeric(nrow(deschutes.m[[2]]))
n0[which(deschutes.m[[2]]$branch == 7)][1] <- 100

# run
d.drift <- stream.drift(n0, l = 100, r = 2, k = 500, Tf = 50, in.mat = deschutes.m$matrix, in.xs = deschutes.m$xs$xs)

#plot
plot <- plot.stream.drift(drift.dat = d.drift, edge.dat = deschutes, xs = deschutes.m$xs)

# sample inds
s.brchs <- c(9, 5)
s.dists <- c(3, 15)
sites <- cbind(s.brchs, s.dists)
ns <- c(100, 100)
inds <- DriftToInds(d.drift, deschutes.m$xs, sites, ns)


# get 2D afs



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
