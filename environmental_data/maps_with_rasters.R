###########################################
############### MAKE MAPS #################
###########################################

rm(list=ls())
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library("rcartocolor")
library("ggthemes")
library(tidyverse)
library(gridExtra)
library(scales)
library(giscoR)

load("rasters4maps.Rdata")

coast <- gisco_get_coastallines(resolution = 3)
load("latlon.Rdata")
latlon = latlon %>% distinct()
latlon$site = substr(rownames(latlon),1,2)

Group1<-c("FW","MC","JT","AH","PT","KR","FB")
Group2<-c("DB","SH","PP","HI","ED","SM","EB")
Group3<-c("LP","PV","SB","DI","SQ")

latlon$clust = ifelse(latlon$site %in% Group1, "SJF",
                      ifelse(latlon$site %in% Group2,"WB","SPS"))



#kd490 sp --------------
kd490_spdf <- as(kd490.r, "SpatialPixelsDataFrame")
kd490_spdf <- as.data.frame(kd490_spdf)
colnames(kd490_spdf) <- c("value", "x", "y")

ggplot() +  
  geom_tile(data=kd490_spdf, aes(x=x, y=y, fill=value)) + 
  scale_fill_viridis_c(option = "mako",direction = -1)+
  coord_equal() +
  theme_map() +
  geom_sf(data=coast,lwd=.3,fill="gray") +
  coord_sf(
    xlim = c(-123.8, -122),
    ylim = c(47.1, 48.25)
  ) +
  #geom_point(data=latlon,aes(lon,lat),color="white",size=4.5)+
  geom_point(data=latlon,aes(lon,lat),color="black",size=6)+
  geom_point(data=latlon,aes(lon,lat,color=clust),size=4.5)+
  scale_color_manual(values = c("hotpink2","darkorange3","orchid4"))+
  theme_classic()+
  xlab("")+
  ylab("")

#ggsave("figures/all5_vars/kd490.png", width = 15, height = 15, units = "cm")
#ggsave("figures/all5_vars/kd490.pdf", width = 15, height = 15, units = "cm")



#kd490w ----------------
library(scales)
kd490winter_spdf <- as(kd490winter.r, "SpatialPixelsDataFrame")
kd490winter_spdf <- as.data.frame(kd490winter_spdf)
colnames(kd490winter_spdf) <- c("value", "x", "y")

ggplot() +  
  geom_tile(data=kd490winter_spdf, aes(x=x, y=y, fill=value)) + 
  scale_fill_viridis_c(option = "mako",direction = -1,limits=c(min(kd490winter_spdf$value),max(kd490_spdf$value)),oob=scales::squish)+
  coord_equal() +
  theme_map() +
  #theme(legend.position="none")+
  geom_sf(data=coast,lwd=.3,fill="gray") +
  coord_sf(
    xlim = c(-123.8, -122),
    ylim = c(47.1, 48.25)
  ) +
  geom_point(data=latlon,aes(lon,lat),color="black",size=6)+
  geom_point(data=latlon,aes(lon,lat,color=clust),size=4.5)+
  scale_color_manual(values = c("hotpink2","darkorange3","orchid4"))+
  theme_classic()+
  xlab("")+
  ylab("")

#ggsave("figures/all5_vars/kd490winter.png", width = 15, height = 15, units = "cm")
#ggsave("figures/all5_vars/kd490winter.pdf", width = 15, height = 15, units = "cm")



#SST W ----
load("sst.wint.raster.Rdata")
sst.wint.df = data.frame("x" = sst.wint.df$x,"y" = sst.wint.df$y,"z" = sst.wint.df$value)
sst.wint.df = sst.wint.df[sst.wint.df$y <= 48.25,]
sst.wint.r = rasterFromXYZ(sst.wint.df)

sst.wint_spdf <- as(sst.wint.r, "SpatialPixelsDataFrame")
sst.wint_spdf <- as.data.frame(sst.wint_spdf)
colnames(sst.wint_spdf) <- c("value", "x", "y")


ggplot() +  
  geom_tile(data=sst.wint_spdf, aes(x=x, y=y, fill=value)) + 
  #scale_fill_gradientn(colours = sst.wint.colours[c(1, seq_along(sst.wint.colours), length(sst.wint.colours))])+
  #scale_fill_viridis_c(direction=1)+
  #paletteer::scale_fill_paletteer_c("grDevices::Purple-Yellow")+
  scale_fill_carto_c(palette = "BurgYl",direction=1)+
  coord_equal() +
  theme_map() +
  geom_sf(data=coast,lwd=.3,fill="gray") +
  coord_sf(
    xlim = c(-123.8, -122),
    ylim = c(47.1, 48.25)
  ) +
  #geom_point(data=latlon,aes(lon,lat),color="white",size=4.5)+
  geom_point(data=latlon,aes(lon,lat),color="black",size=6)+
  geom_point(data=latlon,aes(lon,lat,color=clust),size=4.5)+
  scale_color_manual(values = c("hotpink2","darkorange3","orchid4"))+
  theme_classic()+
  xlab("")+
  ylab("")


#ggsave("figures/all5_vars/sstwinter.png", width = 15, height = 15, units = "cm")

#salinity ----
salinity_spdf <- as(sal.r, "SpatialPixelsDataFrame")
salinity_spdf <- as.data.frame(salinity_spdf)
colnames(salinity_spdf) <- c("value", "x", "y")

min(salinity_spdf$value)
hist(salinity_spdf$value)

salinity_spdf[salinity_spdf$value <= 20,] <- 20

salinity_spdf[which(salinity_spdf$value <= 20),]$value = rescale(salinity_spdf[which(salinity_spdf$value <= 20),]$value, to=c(15,20))


ggplot() +  
  geom_tile(data=salinity_spdf, aes(x=x, y=y, fill=value)) + 
  #scale_fill_gradientn(colours = sal.colours[c(1, seq_along(sal.colours), length(sal.colours))])+
  scale_fill_viridis_c(direction = -1)+
  coord_equal() +
  theme_map() +
  theme(legend.position="none")+
  geom_sf(data=coast,lwd=.3,fill="gray") +
  coord_sf(
    xlim = c(-123.8, -122),
    ylim = c(47.1, 48.25)
  ) +
  #geom_point(data=latlon,aes(lon,lat),color="white",size=4.5)+
  geom_point(data=latlon,aes(lon,lat),color="black",size=6)+
  geom_point(data=latlon,aes(lon,lat,color=clust),size=4.5)+
  scale_color_manual(values = c("hotpink2","darkorange3","orchid4"))+
  theme_classic()+
  xlab("")+
  ylab("")


#ggsave("figures/all5_vars/salinity.png", width = 15, height = 15, units = "cm")


#par ----
par_spdf <- as(par.r, "SpatialPixelsDataFrame")
par_spdf <- as.data.frame(par_spdf)
colnames(par_spdf) <- c("value", "x", "y")


ggplot() +  
  geom_tile(data=par_spdf, aes(x=x, y=y, fill=value)) + 
  #scale_fill_gradientn(colours = par.colours[c(1, seq_along(par.colours), length(par.colours))])+
  scale_fill_viridis_c(option = "magma")+
  coord_equal() +
  theme_map() +
  theme(legend.position="none")+
  geom_sf(data=coast,lwd=.3,fill="gray") +
  coord_sf(
    xlim = c(-123.8, -122),
    ylim = c(47.1, 48.25)
  ) +
  #geom_point(data=latlon,aes(lon,lat),color="white",size=4.5)+
  geom_point(data=latlon,aes(lon,lat),color="black",size=6)+
  geom_point(data=latlon,aes(lon,lat,color=clust),size=4.5)+
  scale_color_manual(values = c("hotpink2","darkorange3","orchid4"))+
  theme_classic()+
  xlab("")+
  ylab("")


#ggsave("figures/all5_vars/par.png", width = 15, height = 15, units = "cm")




