library(sp)
library(rgdal)
library(raster)
library(RColorBrewer)
library(rgeos)
library(maptools)

#stream and watershed infiles
rivers <- shapefile("Project/rivers/rivers.shp")
ws <- shapefile("~/../Downloads/Oregon_Watershed_Councils_2014/Oregon_Watershed_Councils_2014.shp")
dcws <- ws[ws$altName == "Upper Deschutes WC" | ws$altName == "Crooked River WC",]
#plot(dcws)
dcws <- aggregate(dcws)
newCRS <- CRS("+proj=lcc +lat_1=43 +lat_2=45.5 +lat_0=41.75 +lon_0=-120.5 +x_0=399999.9999984001 +y_0=0 +datum=NAD83 +units=ft +no_defs +ellps=GRS80 +towgs84=0,0,0")

#elevation infiles, crop them, fix their row ids, and rbind them.
e_c <- shapefile("Project/Contour_100ce/Contour_100ce.shp")
e_s <- shapefile("Project/Contour_100s/Contour_100s.shp")
e_w <- shapefile("Project/Contour_100w/Contour_100w.shp")
e_c_c <- crop(e_c, dcws)
e_s_c <- crop(e_s, dcws)
e_w_c <- crop(e_w, dcws)
e_c_c <- spChFIDs(e_c_c, as.character(1:nrow(e_c_c)))
at.row <- nrow(e_c_c)
e_s_c <- spChFIDs(e_s_c, as.character((1+at.row):(at.row + nrow(e_s_c))))
at.row <- at.row + nrow(e_s_c)
e_w_c <- spChFIDs(e_w_c, as.character((1+at.row):(at.row + nrow(e_w_c))))
elev <- spRbind(e_c_c, e_s_c)
elev <- spRbind(elev, e_w_c)
elev <- spTransform(elev, newCRS)
remove(e_c, e_s, e_w, e_c_c, e_s_c, e_w_c)


#crop rivers
rivers <- spTransform(rivers, newCRS)
dr <- crop(rivers, dcws)
remove(rivers)

plot(temp, col = "red")
plot(dcws, add = TRUE)
plot(soi, add = TRUE)
plot(elev[as.numeric(elev$CONTOUR) == 2000 | as.numeric(elev$CONTOUR) == 2100,], add = TRUE, col = "grey")

#plot(dr)
#plot(dcws, add = TRUE)

soi <- dr[dr$NAME == "Pilot Butte Canal" |
            dr$NAME == "Deschutes River" |
            dr$NAME == "Crooked River"|
            dr$STREAMS_ID == "18191",]
soi <- disaggregate(soi)

#plot(soi, col = "red")
#plot(dcws, add = TRUE)
#plot(elev[as.numeric(elev$CONTOUR) %% 500 == 0,], col = brewer.pal(9, "Greys"), add = TRUE)

test <- GIS.to.Edge(soi, ei = elev, ei.c = "CONTOUR", TRUE)

test$branch.length <- test$branch.length/5200
test$start.weight <- 1
test$end.weight <- 1
test$branch.length <- round(test$branch.length)

test.m <- Edge.to.Matrix(test, 2, 5, 1)
#normalize
for (i in 1:ncol(test.m[[1]])){
  test.m[[1]][,i] <- test.m[[1]][,i]/sum(test.m[[1]][,i])
}


#barriers <- shapefile("Project/ODFW_44_5_ofpbds_shp/ofpbds_pt.shp")
#barriers <- spTransform(barriers, newCRS)
#barriers <- crop(barriers, dcws)
#impass <- barriers[barriers$fpbFPasSta == "Blocked"|barriers$fpbFPasSta == "Partial",]
#plot(impass)
#plot(soi, add = TRUE)

#wbs <- shapefile("Project/waterbodies/waterbodies.shp")
#wbs <- spTransform(wbs, newCRS)
#plot(soi)
#plot(impass, add = TRUE)
#plot(wbs, add = TRUE)

#sn <- SpatialLinesNetwork(soi, 20000)

#plot(sn@g)
#plot(sn@g$x, sn@g$y, col = sn@g$n, pch = 16, cex = 2, asp = 1)
#lines(soi)
#text(sn@g$x, sn@g$y, E(sn@g), pos = 4)





rownames(zd_es) <- 1:nrow(zd_es)
temp <- SpatialPoints(zd_es, proj4string = newCRS)
v.els <- numeric(nrow(data.frame(temp)))
for (i in 1:length(shortest.dists)) {
  print(i)
  v.els[i] <- as.numeric(elev$CONTOUR[which.min(gDistance(temp[i,], elev, byid = TRUE))])
}
zd_es <- cbind(zd_es, v.els)

as.numeric(elev$CONTOUR[which.min(gDistance(temp[1,], elev, byid = TRUE))])

ds <- gDistance(temp[1,], elev, byid = TRUE)
ds <- cbind(1:length(ds), ds)
colnames(ds) <- c("index", "dist")
ds <- as.data.frame(ds)
ds <- arrange(ds, dist)
ds <- ds[!duplicated(ds[,2]),]
w1 <- 1-(ds[1,2]/(ds[1,2]+ds[2,2]))
w2 <- 1-w1
wd <- w1*as.numeric(elev$CONTOUR[ds[1,1]]) + w2*as.numeric(elev$CONTOUR[ds[2,1]])
