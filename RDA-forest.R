

#Nereo.RDAforest.Rdata

rm(list=ls())
# edit this line to point to the path you downloaded the files to:

library(RDAforest)
library(ggplot2)
library(gradientForest)
library(mapdata)
library(vegan)
library(tidyverse)
library(randomForest)
library(giscoR)
library(cluster)
library(clv)
library(sf)
library(rnaturalearth)
library(sp)

source("RDAforest_functions.R") #for GCD

coast <- gisco_get_coastallines(resolution = 3)
landscape=st_as_sf(ne_countries(scale="medium",continent="north america"))
mapdata= ne_countries(scale="large",continent="north america")
theme_set(theme_bw())

# the next four are only to add outline of coasts to the last plot, at the very last line of this script.
# Skip if having trouble installing these packages.
# library(maps)
# library(mapdata)

# loading all data: 
load("RDAforest.input65.Rdata")

# load("merged_rasters_new_sst_ssm.Rdata") #adjusted wint sst
# envc = merge.dist
# envc = envc[,c(2:7,1)]
# envc = envc[seq(1,nrow(envc), 25),]
# colnames(envc) = c("x","y","SalinityAVG2","SST.month.Winter","kd490.Winter","PAR.Summer2","kd490.Spring")
# 

#OLD
#load("merged_rasters_sst_fill.Rdata")
load("merged_rast_fresh_kd490.Rdata")
envc=merge.dist
envc = envc[,c(2:7,1)]
envc = envc[seq(1,nrow(envc), 25),]
colnames(envc) = c("x","y","kd490.Winter","PAR.Summer2","SST.month.Winter","SalinityAVG2","kd490.Spring")

#trim 
env = env[,-c(4,7,8,9)]

plot(ord$CA$eig/sum(ord$CA$eig),xlab="PC",ylab="proportion of variance explained")

so=data.frame(scores(ord,scaling=1,display="sites"))
so$ecotype = ecotype$ecotype
ggplot(so,aes(MDS1,MDS2,color=ecotype))+geom_point()+coord_equal()+theme_classic()+scale_color_manual(values = c("red","blue","green"))

GCD=gcd.dist(latlon)
latlon.gcd=GCD[[1]]
distGCD=GCD[[2]]

# water.dist = read.table("pairwise_over_water_distances.txt",header=T)
# water.dist.mat=xtabs(Distance~.,water.dist)

plot(as.dist(cordist)~distGCD,pch=16,cex=0.6,col=rgb(0,0,0,alpha=0.2))

protest(capscale(distGCD~1),capscale(cordist~1))

# getting the first two PCs of log-distance matrix
d.ord=capscale(distGCD~1)
pcs.d=scores(d.ord,scaling=1,display="sites",choices=c(1:2))

# constructing and plotting partial ordination
ord1=capscale(cordist~1+Condition(as.matrix(pcs.d)))
plot(ord1$CA$eig/sum(ord1$CA$eig),xlab="PC",ylab="proportion of variance explained")

so1=data.frame(scores(ord1,scaling=1,display="sites"))
ggplot(so1,aes(MDS1,MDS2,color=ecotype$ecotype))+geom_point()+coord_equal()+theme_bw()


pc=hclust(as.dist(1-cor(env)))
plot(pc)
abline(h=0.1, col="red")

#keep PCs
gf=makeGF(ord1,env,pcs2keep=c(1:20))
# predicted PCs, and how well they are predicted (per-PC R2s)
gf$result
# MDS1        MDS2        MDS3        MDS4        MDS5        MDS6        MDS7       MDS12       MDS13       MDS19 
# 0.926802932 0.912709139 0.762349183 0.655435075 0.622464907 0.276246431 0.764447929 0.040492875 0.493300159 0.00822705

#how much variance does our model capture?
# rescaling to proportion of total variance (based on eigenvalues in ord1)
eigen.var=(ord1$CA$eig/sum(ord1$CA$eig))[names(gf$result)]
# total variance explained by model
sum(eigen.var*gf$result)
#[1] 0.2348012

# setting the number of PCs to keep
tokeep=5
# computing properly scaled importances:
imps=data.frame(importance_RDAforest(gf,ord1))
# some data frame housekeeping...
names(imps)="R2"
imps$var=row.names(imps)
# reordering predictors by their importances:
imps$var=factor(imps$var,levels=imps$var[order(imps$R2)])
# plotting
ggplot(imps,aes(var,R2))+geom_bar(stat="identity")+coord_flip()+theme_bw()

#Let’s plot importances of all predictors. These are properly scaled to correspond to the proportion of variance 
#explained by the predictor in the whole dataset.

#turnover curves
gf3=makeGF(ord1,env,pcs2keep=c(1:3))
plot_gf_turnovers(gf3,imps$var[1:5])

cordist = as.data.frame(cordist)
rownames(cordist) = rownames(latlon)
colnames(cordist) = rownames(latlon)

mm=mtrySelJack(Y=cordist,X=env,covariates=pcs.d,nreps=50,prop.positive.cutoff=0.3, top.pcs=tokeep)
mm$goodvars

ggplot(mm$prop.positive,aes(var,prop.positive))+
  geom_bar(stat="identity")+
  coord_flip()+
  geom_hline(yintercept=0.3,col="red")

ggplot(mm$delta,aes(var,values))+
  geom_boxplot(outlier.shape = NA)+
  coord_flip()+
  geom_hline(yintercept=0,col="red")


#try getting rid of nearshore weirdness
#trims = rownames(envc[envc$x >= -122.55 & envc$y >= 48.1,])
#merge.trim = envc[-which(rownames(envc) %in% trims),]

oj=ordinationJackknife(Y=cordist,X=env[,mm$goodvars],newX=envc,covariates=pcs.d,nreps=50,top.pcs=tokeep,extra=0.1)

sum(oj$median.importance)

par(mfrow = c(2,3))
for(i in 1:length(mm$goodvars)){
  plot_turnover(oj,envc[oj$goodrows,],names(oj$median.importance)[i])
}

dev.off()

ggplot(oj$all.importances,aes(variable,importance))+geom_boxplot(outlier.shape = NA)+coord_flip()

envcbars = oj$all.importances
lsbars = envcbars %>% group_split(variable)
ls2means = list()
for (i in seq_along(lsbars)) {
  ls2means[[i]] = data.frame("var" = unique(lsbars[[i]]$variable),
                             "mean.imp" = mean(lsbars[[i]]$importance))
}

impbars = bind_rows(ls2means)
ggplot(impbars,aes(var,mean.imp))+geom_bar(stat="identity")+coord_flip()+theme_bw()



#get min/max importance for each var
# sal = oj$all.importances[which(oj$all.importances$variable == "SalinityAVG2"),]
# wint = oj$all.importances[which(oj$all.importances$variable == "SST.month.Winter"),]
# 
# min(sal$importance) #0.06
# max(sal$importance) #0.08
# 
# min(wint$importance) #0.03
# max(wint$importance) #0.05

#ADAPTIVE NEIGHBORHOODS -------------
#Forming predictions based on averaging replicates from ordinationJackknife
# spots on the map that are within modeled parameter range:
goods=oj$goodrows

# predictor data restricted to only those spots:
ras2=envc[which(goods),]
xy2=envc[goods,c("x","y")]
names(xy2)=c("lon","lat")
rfpreds=oj$predictions.direct
turnovers=oj$predictions.turnover
bests=names(oj$median.importance)

plot_adaptation(rfpreds,ras2[,bests],xy2,main="unclustered",
                # options affecting PCA plot:
                rangeExp=1.5,
                scal=20,
                jitscale=0.05,
                # options affecting colors:
                color.scheme="011",
                lighten=0.8
)

plot(coast$geometry, add = TRUE, col = "grey")
points(x=latlon$lon,y=latlon$lat)

#clustering adaptive neighborhoods -----------
#GET THE BOUNDARIES OF THESE NEIGHBORHOODS
#
#Now our goal is to break the continuous colors in the map above into “adaptive neighborhoods” - 
#bounded areas likely to contain similarly-adapted organisms

load("rdaForest_weird_blue.RData")

#method 1: cluster spatial points based on random forest predictions
pa1=plot_adaptation(rfpreds,ras2[,bests],xy2,main="direct preds",
                    # options affecting PCA plot:
                    rangeExp=1,
                    scal=30,
                    jitscale=0.05,
                    # options affecting map and PCA colors:
                    color.scheme="011",
                    lighten=0.8,
                    # options affecting clustering:
                    cluster.guide = NULL,
                    nclust=20,
                    cluster.merge.level=0.35
)


pa2=plot_adaptation(rfpreds,ras2[,bests],xy2,main="turnovers",
                    # options affecting PCA plot:
                    rangeExp=1.8,
                    scal=15,
                    jitscale=2,
                    # options affecting map and PCA colors:
                    color.scheme="011",
                    lighten=0.8,
                    # options affecting clustering:
                    cluster.guide = turnovers,
                    nclust=20,
                    cluster.merge.level=0.35
)

plot(coast$geometry, add = TRUE, col = "grey")

#nice map (clustered)
plot_nice_map(xy2,mapdata=landscape,map.on.top=T,size=1,cols=pa2$colors,
              overlay.points = latlon,
              cols.points = ecotype$color,size.points=2)


ecotype$color = ifelse(ecotype$ecotype == "SJF", "hotpink2",
                       ifelse(ecotype$ecotype == "WB","orchid4","darkorange3"))
colors2plot = cbind(latlon,ecotype)
colors2plot = colors2plot[,c(1,2,5,6)]
colors2plot = colors2plot %>% distinct()

# plot(x = xy2$lon,y=xy2$lat,col=pa2$colors,xlab=c("longitude"),ylab=c("latitude"),pch=16,xlim=c(-124,-122.3))
plot(x = xy2$lon,y=xy2$lat,col=pa2$colors,xlab=c("longitude"),ylab=c("latitude"),pch=16,xlim=c(-123.8,-122.3)) #new xlim
plot(coast$geometry, add = TRUE, col = "lightgrey")
points(x=colors2plot$lon,y=colors2plot$lat,col = "black",pch=1,cex=1.6)
points(x=colors2plot$lon,y=colors2plot$lat,col = "black",pch=1,cex=1.7)
points(x=colors2plot$lon,y=colors2plot$lat,col = "black",pch=1,cex=1.8)
points(x=colors2plot$lon,y=colors2plot$lat,col = alpha("white",0.6),pch=16,cex=1.5)
points(x=colors2plot$lon,y=colors2plot$lat,col = alpha(colors2plot$color,0.8),pch=16,cex=1.5)

#save(oj,mm,tokeep,xy2,gf,pa2,file=paste("rdaForest_result.RData"))
#save(oj,mm,tokeep,xy2,gf,pa2,file=paste("rdaForest_result_ssm.RData"))
#save(oj,mm,tokeep,xy2,gf,pa2,file=paste("rdaForest_weird_blue.RData"))


#env mismatch ----------
rm(list=ls())
#load data
source("RDAforest_functions.R")
coast <- gisco_get_coastallines(resolution = 3)
load("RDAforest.input65.Rdata")
load("merged_rasters_sst_fill.Rdata")
envc = merge.dist
envc = envc[,c(2:9,1)]
envc = envc[seq(1,nrow(envc), 25),]
colnames(envc) = c("x","y","SST.August","SST.month.Winter","SalinityAVG2","PAR.Winter2","kd490.Winter","PAR.Summer2","kd490.Spring")
env = env[,colnames(env) %in% colnames(envc)]
which(colnames(env) %in% colnames(envc))
env = env[,-c(4,7)]
load("rdaForest_result.RData")

#make plots
important=names(oj$median.importance)
ggplot(oj$all.importances[oj$all.importances$variable %in% important,],aes(variable,importance))+
  geom_boxplot(outlier.shape = NA)+
  coord_flip()+
  theme_bw()

sum(oj$median.importance)

#get SPS coords
sps.ll = data.frame(latlon,"site" = substr(rownames(latlon),1,2))
sps.ll = sps.ll[sps.ll$site == "DI" | sps.ll$site == "LP" | sps.ll$site == "PV" | sps.ll$site == "SB" | sps.ll$site == "SQ",]
sps.ll = sps.ll %>% distinct()
plot(sps.ll$lon,sps.ll$lat)
text(sps.ll$site,x=sps.ll$lon,y=sps.ll$lat)

plot(latlon$lon,latlon$lat)
text(substr(rownames(latlon),1,2),x=latlon$lon,y=latlon$lat)

# reef=c(-122.363,48.051)
# reef=c(-122.901,47.169) #SQ
# reef=c(-122.565,47.245) #DI
# reef=c(-122.545,47.949) #DB
# reef=c(-122.431,47.681) #SM
# reef=c(-122.535,47.303) #SB
reef=c(-123.051,48.123) #JT

rows2keep = names(which(oj$goodrows == TRUE))
XY = envc[which(rownames(envc) %in% rows2keep),]
xy.ll = XY[,1:2]
XY = XY[,which(colnames(XY) %in% mm$goodvars)]
XY = cbind(xy.ll,XY)

library(terra)
library(viridis)
xy2 = XY
rf.preds=terra::extract(lonlat2raster(xy2,oj$predictions.direct),data.frame(rbind(reef)),method="simple")[,-1]
sc=adapt_scale(oj$predictions.direct)[2]
# computing distances between future-needed and present-day adaptation
agf=env_mismatch(X=rf.preds,Y=oj,sy=XY,sc=sc)

# agf1=as.data.frame(agf)
# which(is.na(rowSums(agf1)))

ggplot(agf,aes(x,y,color=env.mismatch))+
  geom_point(size=1,alpha=0.5)+
  #geom_point(size = 2.5,alpha = .5, stroke=NA, position=position_jitter(height=.009, width=.009))+
  coord_equal()+
  theme_minimal()+
  #ylim()+
  #xlim(-125,-121)+
  scale_color_viridis()+
  geom_point(data=data.frame(rbind(reef)),aes(X1,X2),pch=18,size=5,col="red")

agf.r = agf[,c(1,2,8)]
agf.r = rasterFromXYZ(agf.r)

r = agf.r
plot(r)
wm<-focalWeight(r,d=res(r)[1]*1.001,type="circle", fillNA = TRUE) # original weight matrix
w<-wm
w[w>0]<-1
# r1<-focal(r,w=w,fun=mean,na.rm=T)
r1<-raster::focal(r,w=w,fun=mean,na.rm=T,NAonly=T)
plot(r1,main="focal weight + mean")
#plot(coast$geometry, add = TRUE, col = "grey")

agf_spdf <- as(r1, "SpatialPixelsDataFrame")
agf_spdf <- as.data.frame(agf_spdf)
colnames(agf_spdf) <- c("value", "x", "y")

latlon.sm = data.frame("lat" = latlon$lat,"lon" = latlon$lon,"site" = substr(rownames(latlon),1,2)) %>% distinct()

#ll2plot = latlon.sm[-c(16),]

ggplot()+
  geom_raster(data=agf_spdf, aes(x=x, y=y, fill=value),interpolate = T)+
  scale_fill_viridis()+
  coord_equal() +
  theme_map() +
  geom_sf(data=coast,lwd=.3,fill="gray") +
  coord_sf(
    xlim = c(-123.8, -122),
    ylim = c(47.1, 48.25)
  ) +
  #theme(legend.position = "none")+
  geom_point(data=latlon,aes(lon,lat),color="black",size=4,alpha=0.9)+
  geom_point(data=data.frame(rbind(reef)),aes(X1,X2),pch=18,size=5,col="red")
#geom_text(data=data.frame(rbind(reef)),aes(X1,X2),label="★", size=9, family = "HiraKakuPro-W3",color="black")+
#geom_text(data=data.frame(rbind(reef)),aes(X1,X2),label="★", size=6, family = "HiraKakuPro-W3",color="yellow")

ggplot()+
  geom_tile(data=ras_df, aes(x=x, y=y, fill=value))+
  scale_fill_viridis()+
  coord_equal() +
  theme_map() +
  geom_sf(data=coast,lwd=.3,fill="gray") +
  coord_sf(
    xlim = c(-123.8, -122),
    ylim = c(47.1, 48.25)
  ) +
  #theme(legend.position = "none")+
  geom_point(data=latlon,aes(lon,lat),color="black",size=4,alpha=0.9)+
  geom_point(data=data.frame(rbind(reef)),aes(X1,X2),pch=18,size=5,col="red")


ras_df <- agf_spdf %>% mutate(across(c(x, y), round, digits = 2))

ggplot()+
  geom_raster(data=ras_df, aes(x=x, y=y, fill=value),interpolate = T)+
  scale_fill_viridis()+
  coord_equal() +
  theme_map() +
  geom_sf(data=coast,lwd=.3,fill="gray") +
  coord_sf(
    xlim = c(-123.8, -122),
    ylim = c(47.1, 48.25)
  ) +
  #theme(legend.position = "none")+
  geom_point(data=latlon,aes(lon,lat),color="black",size=4,alpha=0.9)+
  geom_point(data=data.frame(rbind(reef)),aes(X1,X2),pch=18,size=5,col="red")
#geom_text(data=data.frame(rbind(reef)),aes(X1,X2),label="★", size=9, family = "HiraKakuPro-W3",color="black")+
#geom_text(data=data.frame(rbind(reef)),aes(X1,X2),label="★", size=6, family = "HiraKakuPro-W3",color="yellow")


agf.round <- agf %>% mutate(across(c(x, y), round, digits = 2))

ggplot()+
  geom_raster(data=agf.round, aes(x=x, y=y, fill=env.mismatch),interpolate = T)+
  scale_fill_viridis()+
  coord_equal() +
  theme_map()+
  geom_sf(data=coast,lwd=.3,fill="gray") +
  coord_sf(
    xlim = c(-123.8, -122),
    ylim = c(47.1, 48.25)
  ) +
  geom_point(data=data.frame(rbind(reef)),aes(X1,X2),pch=18,size=5,col="red")



