
#create node list to use as imput for salinity_ssm.py

library(ncdf4)
library(CFtime)
library(lattice)

#a closer look at that bathy
bathy = read.delim("SSM_Original_Bathy_UTM_Zone10_NAVD88.xyz",sep = " ")
bathy$node = as.numeric(rownames(bathy))
psnodes = read.table("~/Downloads/PugetNodeList.csv")
psnodes$V1 = as.numeric(psnodes$V1)

ncin <- nc_open("HYD_SSM_2014_00301.nc")
ncsmall = ncin$var$X$dim[[1]]

length(which(ncsmall$vals %in% psnodes$V1))
bathy$col = ifelse(rownames(bathy) %in% psnodes$V1,"red","black")
pscols = bathy[bathy$col == "red",]

plot(bathy$XYZ,bathy$X.16012)
points(pscols$XYZ,pscols$X.16012,col="red")

#what are the z values?
min(pscols$points.)
max(pscols$points.)

min(bathy$points.)
max(bathy$points.)

hist(bathy$points.,breaks=1000)
hist(pscols$points.,breaks=1000)

length(which(bathy$node %in% pscols$node))

#shallow = bathy[bathy$points. <= 284.052,]
shallow = bathy
plot(shallow$XYZ,shallow$X.16012)
length(which(pscols$node %in% shallow$node))

max(shallow$XYZ)
min(shallow$XYZ)

max(pscols$XYZ)
min(pscols$XYZ)

#remove weird lower bits
plot(pscols$XYZ,pscols$X.16012)
abline(h=5150000,col="red")
pscols = pscols[pscols$X.16012 >= 5150000,]
plot(pscols$XYZ,pscols$X.16012)
length(which(pscols$node %in% shallow$node))


#extend to sjf
plot(shallow$XYZ,shallow$X.16012)
points(pscols$XYZ,pscols$X.16012,col="red")
abline(h=5380000,col="red")
abline(v=380000,col="red")
abline(h=min(pscols$X.16012))

region = shallow[shallow$XYZ > 380000 &
                   shallow$X.16012 < 5380000 &
                   shallow$X.16012 > min(pscols$X.16012) ,]

plot(region$XYZ,region$X.16012)
abline(v=425000, col="red")
abline(h=5310000)

region = region[-which(region$XYZ < 420000 & region$X.16012 < 5310000),]

plot(region$XYZ,region$X.16012)
points(pscols$XYZ,pscols$X.16012,col="red")

length(which(pscols$node %in% region$node))

#make sure all nodes are in original bathy
length(which(region$node %in% bathy$node))
`%nin%` <- Negate(`%in%`)
which(region$node %nin% bathy$node)

nodes = as.numeric(region$node)

max(nodes)
min(nodes)

max(bathy$node)
node = nodes[1:7199] #for some reason the last node isn't in the nc data for some files

#save(region,file="nodesinregion.Rdata")
#write.table(nodes,file="PugetNodeList.csv",sep = ",",col.names = F,row.names = F,quote = F)

#now figure out how to make a raster
library(RColorBrewer)
library(tidyverse)

rm(list=ls())
load("nodesinregion.Rdata")
result=read.table("avgMinDailySal.csv",header = T,sep=",")
colnames(result) = c("node","val")

all = left_join(region,result,by="node")


cols = brewer.pal(4, "Blues")
pal = colorRampPalette(cols)

#plot(all$XYZ,all$X.16012,col=pal(nrow(all)[all$val]))

plot(all$XYZ,all$X.16012,col=all$val)

ggplot(all,aes(XYZ,X.16012,color=val))+
  geom_point(alpha=0.5)+
  #scale_color_continuous()
  scale_colour_gradient(low = "red", high = "black")+
  theme_minimal()

res = data.frame("x" = all$XYZ,"y"=all$X.16012,"z" = all$val)
#res = na.omit(res)
write.table(res,file="salinity_results.csv",col.names = T,row.names = F,quote = F,sep = ",")

#interpolate in arcgis

#import from arcgis ----------
Salinity<-stack("sal_arcgis/sal7.tif")
r = Salinity[[1]]
plot(Salinity)
#plot(coast$geometry, add = TRUE, col = "grey")
#points(latlon$lon,latlon$lat)

#make sure we've got data at these points
coords = latlon[,c(2,1)]
colnames(coords) = c("Longitude","Latitude")
coordinates(coords) <-  ~Longitude+Latitude

x_coords = c(-122,-122,-124, -124, -122.7,  -122.65, -122.6, -123.4,-122)
y_coords = c(46.9,48.5,48.5, 47.9, 47.9,  47.6 ,47.5,46.9,46.9)

poly1 <- sp::Polygon(cbind(x_coords,y_coords)) #make polygon from vertices matrix
firstPoly <- sp::Polygons(list(poly1), ID = "A") #make polygon class
str(firstPoly,1)
firstSpatialPoly <- sp::SpatialPolygons(list(firstPoly))
plot(poly1@coords)
plot(firstSpatialPoly)

valscoords = as.data.frame(raster::extract(r,firstSpatialPoly,cellnumbers=T))
sr <- cbind(valscoords, xyFromCell(r,valscoords[,1]))
sr$value[sr$value <= 0] <- NA

sr = rasterFromXYZ(sr[,c(3,4,2)])
plot(sr)
plot(coast$geometry, add = TRUE, col = "grey")


wm<-focalWeight(r,d=res(r)[1]*1.001,type="circle", fillNA = TRUE) # original weight matrix
w<-wm
w[w>0]<-1
r1<-focal(sr,w=w,fun=mean,na.rm=T)
plot(r1,main="focal weight + mean")
#plot(coast$geometry, add = TRUE, col = "grey")



#increase res
diss <- disaggregate(r1, fact=25,method="bilinear")
plot(diss,main=paste("fact"))
plot(coast$geometry, add = TRUE, col = "grey")

valscoords = as.data.frame(raster::extract(diss,firstSpatialPoly,cellnumbers=T))
sal.points = as.data.frame(raster::extract(r1,coords,cellnumbers=T))
sal.avg.df <- cbind(valscoords, xyFromCell(diss,valscoords[,1]))
save(sal.points,sal.avg.df,file="sal.raster.Rdata")
