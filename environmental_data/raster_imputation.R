
#goal: impute missing data and standardize cell size across all rasters
#note: salinity is imputed separately in arcGIS because it is sourced from the salish sea model

library(raster)
library(sp)
library(giscoR)

rm(list=ls())
load("Salinity_rasters.Rdata")
load("PAR_rasters.Rdata")
load("kd490_rasters.Rdata")

load("env_65.Rdata") #buoy data and site coordinates

coast <- gisco_get_coastallines(resolution = 3) #for plotting nicely

#polygon barrier so we don't impute over land
x_coords = c(-122,-122,-124, -124, -122.7,  -122.65, -122.6, -123.4,-122)
y_coords = c(46.9,48.5,48.5, 47.9, 47.9,  47.6 ,47.5,46.9,46.9)

poly1 <- sp::Polygon(cbind(x_coords,y_coords)) #make polygon from vertices matrix
firstPoly <- sp::Polygons(list(poly1), ID = "A") #make polygon class
str(firstPoly,1)
firstSpatialPoly <- sp::SpatialPolygons(list(firstPoly))
plot(poly1@coords)
plot(firstSpatialPoly)

#get points for more extracts
env = env[rownames(latlon),]
coords = latlon[,c(2,1)]
colnames(coords) = c("Longitude","Latitude")
coordinates(coords) <-  ~Longitude+Latitude

#imputation and resolution
f = 30

# #SST winter ---------------
#impute missing data
r = SST.winter
#try calculating a matrix...?
wm<-focalWeight(r,d=res(r)[1]*1.001,type="circle", fillNA = TRUE) # original weight matrix
w<-wm
w[w>0]<-1
r1<-focal(r,w=w,fun=mean,na.rm=T)
plot(r1,main="focal weight + mean")
plot(coast$geometry, add = TRUE, col = "grey")
#increase res
diss <- disaggregate(r1, fact=f,method="bilinear")
plot(diss,main=paste("fact",f))
plot(coast$geometry, add = TRUE, col = "grey")
#sst.winter.impute = diss

valscoords = as.data.frame(extract(diss,firstSpatialPoly,cellnumbers=T))
sst.wint.df <- cbind(valscoords, xyFromCell(diss,valscoords[,1]))
# 
# plot(x=sst.wint.df$x,y=sst.wint.df$y)
# plot(coast$geometry, add = TRUE, col = "grey")
# points(x=x_coords,y=y_coords,col = c("red"),pch=16)
# 
# 
# max(env$SST.month.Winter) #9.143077
# min(env$SST.month.Winter) #7.726121
# 
# max(na.omit(sst.wint.df$value)) #8.231207
# min(na.omit(sst.wint.df$value)) #7.583944

#par summer  ---------------
r = par.summer
wm<-focalWeight(r,d=res(r)[1]*1.001,type="circle", fillNA = TRUE) # original weight matrix
w<-wm
w[w>0]<-1
r1<-focal(r,w=w,fun=mean,na.rm=T)
plot(r1,main="focal weight + mean")
plot(coast$geometry, add = TRUE, col = "grey")
#increase res
diss <- disaggregate(r1, fact=f,method="bilinear")
plot(diss,main=paste("fact",f))
#plot(coast$geometry, add = TRUE, col = "grey")
points(x=latlon$lon,y=latlon$lat)

valscoords = as.data.frame(extract(diss,firstSpatialPoly,cellnumbers=T))
par.summer.df <- cbind(valscoords, xyFromCell(diss,valscoords[,1]))

valscoords = as.data.frame(extract(diss,coords,cellnumbers=T))
#par.summer.points <- cbind(valscoords, xyFromCell(diss,valscoords[,1]))


# max(env$PAR.Summer2) #48.76606
# min(env$PAR.Summer2) #45.73988
# 
# max(na.omit(par.summer.df$value)) #50.19682
# min(na.omit(par.summer.df$value)) #43.10736

#kd490  winter ---------------
r = Winter.Avg.kd490
wm<-focalWeight(r,d=res(r)[1]*1.001,type="circle", fillNA = TRUE) # original weight matrix
w<-wm
w[w>0]<-1
r1<-focal(r,w=w,fun=mean,na.rm=T)
plot(r1,main="focal weight + mean")
plot(coast$geometry, add = TRUE, col = "grey")
#increase res
diss <- disaggregate(r1, fact=f,method="bilinear")
plot(diss,main=paste("fact",f))
#plot(coast$geometry, add = TRUE, col = "grey")
points(x=latlon$lon,y=latlon$lat)

valscoords = as.data.frame(extract(diss,firstSpatialPoly,cellnumbers=T))
kd490.winter.df <- cbind(valscoords, xyFromCell(diss,valscoords[,1]))

# valscoords = as.data.frame(extract(diss,coords,cellnumbers=T))
# kd490.winter.points <- cbind(valscoords, xyFromCell(diss,valscoords[,1]))

#save(kd490.winter.df,file="kd490_raster.imp.Rdata")

max(env$kd490.Winter) #4.060333
min(env$kd490.Winter) #0.1839297

max(na.omit(par.summer.df$value)) #50.19682
min(na.omit(par.summer.df$value)) #43.10736

#
#kd490  spring ---------------
r = Spring.Avg.kd490
wm<-focalWeight(r,d=res(r)[1]*1.001,type="circle", fillNA = TRUE) # original weight matrix
w<-wm
w[w>0]<-1
r1<-focal(r,w=w,fun=mean,na.rm=T)
plot(r1,main="focal weight + mean")
plot(coast$geometry, add = TRUE, col = "grey")
#increase res
diss <- disaggregate(r1, fact=f,method="bilinear")
plot(diss,main=paste("fact",f))
plot(coast$geometry, add = TRUE, col = "grey")
points(x=latlon$lon,y=latlon$lat)

valscoords = as.data.frame(extract(diss,firstSpatialPoly,cellnumbers=T))
kd490.spring.df <- cbind(valscoords, xyFromCell(diss,valscoords[,1]))
#save(kd490.winter.df,file="kd490_raster.imp.Rdata")
# valscoords = as.data.frame(extract(diss,coords,cellnumbers=T))
# kd490.spring.points <- cbind(valscoords, xyFromCell(diss,valscoords[,1]))

max(env$kd490.Winter) #4.060333
min(env$kd490.Winter) #0.1839297

max(na.omit(par.summer.df$value)) #50.19682
min(na.omit(par.summer.df$value)) #43.10736


#
Pattern1<-grep(".points",names(.GlobalEnv),value=TRUE)
list.allvar<-do.call("list",mget(Pattern1))


####################################-----------------------------------------------
###### BIND RASTERS TOGETHER #######
####################################

#RDAforest requires that all rasters be in one dataframe
#x y coordinates for each raster must match up!

load("salinity.raster.Rdata") #load the salinity raster from arcGIS
#sal.avg.df = data.frame("cell" = NA, "x" = sal.avg.df$x,"y"=sal.avg.df$y,"value" = sal.avg.df$value)

Pattern1<-grep(".df",names(.GlobalEnv),value=TRUE)
list.allvar<-do.call("list",mget(Pattern1))

save(list.allvar,file="list.allvar.Rdata") 

#bind together
rm(list = ls())
load("list.allvar.Rdata")

library(tidyverse)
for (i in names(list.allvar)) {
  list.allvar[[i]][,c(3,4)] = round(list.allvar[[i]][,c(3,4)],3)
  list.allvar[[i]] = list.allvar[[i]][,-c(1)]
}

for (i in names(list.allvar)) {
  colnames(list.allvar[[i]]) = c(names(list.allvar[i]),"x","y")
}

merge=list.allvar %>% reduce(full_join, by = c("x","y"))
merge = na.omit(merge)

merge.dist = merge %>% distinct()

length(unique(merge.dist$x))
length(unique(merge.dist$y))

min(merge.dist$y)
merge.dist[which(round(merge.dist$x,2) >= 47.16),]

save(merge.dist,file="rasters.Rdata") #for RDAforest input

