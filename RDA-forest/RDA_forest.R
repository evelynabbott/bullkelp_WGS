


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

#get coastlines for pretty maps
coast <- gisco_get_coastallines(resolution = 3)
landscape=st_as_sf(ne_countries(scale="medium",continent="north america"))
mapdata= ne_countries(scale="large",continent="north america")
theme_set(theme_bw())


# loading all data: 
load("RDAforest.input65.Rdata")
load("rasters.Rdata") 
envc=merge.dist
envc = envc[,c(2:7,1)]
envc = envc[seq(1,nrow(envc), 25),]
colnames(envc) = c("x","y","kd490.Winter","PAR.Summer2","SST.month.Winter","SalinityAVG2","kd490.Spring")

#trim unimportant - keep all variables first time around
env = env[,-c(4,7,8,9)]

plot(ord$CA$eig/sum(ord$CA$eig),xlab="PC",ylab="proportion of variance explained")

so=data.frame(scores(ord,scaling=1,display="sites"))
so$ecotype = ecotype$ecotype
ggplot(so,aes(MDS1,MDS2,color=ecotype))+geom_point()+coord_equal()+theme_classic()+scale_color_manual(values = c("red","blue","green"))

#
#great circle dist-----------------------------------------
# GCD=gcd.dist(latlon)
# latlon.gcd=GCD[[1]]
# distGCD=GCD[[2]]
# 
# # water.dist = read.table("pairwise_over_water_distances.txt",header=T)
# # water.dist.mat=xtabs(Distance~.,water.dist)
# 
# plot(as.dist(cordist)~distGCD,pch=16,cex=0.6,col=rgb(0,0,0,alpha=0.2))
# 
# protest(capscale(distGCD~1),capscale(cordist~1))
# 
# # getting the first two PCs of log-distance matrix
# d.ord=capscale(distGCD~1)
# pcs.d=scores(d.ord,scaling=1,display="sites",choices=c(1:2))
# 
# # constructing and plotting partial ordination
# ord1=capscale(cordist~1+Condition(as.matrix(pcs.d)))
# plot(ord1$CA$eig/sum(ord1$CA$eig),xlab="PC",ylab="proportion of variance explained")
# 
# so1=data.frame(scores(ord1,scaling=1,display="sites"))
# ggplot(so1,aes(MDS1,MDS2,color=ecotype$ecotype))+geom_point()+coord_equal()+theme_bw()


#pairwise marine dist -------------------------
mdist = read.table("pairwise_over_water_distances.txt",header = T)
inds = rownames(env)

pairs = t(as.data.frame(combn(inds, 2)))
pairs = as.data.frame(pairs)
pairs$Pop.i = substr(pairs$V1,1,2)
pairs$Pop.j = substr(pairs$V2,1,2)

all = full_join(pairs,mdist,by=c("Pop.i","Pop.j"))
all[is.na(all)] <- 0

df = all[,c(1,2,5)]
mardist <- with(df, Distance)
nams <- with(df, unique(c(as.character(V1), as.character(V2))))
attributes(mardist) <- with(df, list(Size = length(nams),
                                     Labels = nams,
                                     Diag = FALSE,
                                     Upper = FALSE,
                                     method = "user"))
class(mardist) <- "dist"


plot(as.dist(cordist)~mardist,pch=16,cex=0.6,col=rgb(0,0,0,alpha=0.2))
protest(capscale(mardist~1),capscale(cordist~1))

#Procrustes Sum of Squares (m12 squared):        0.2917  #aka: dissimilarity
#Correlation in a symmetric Procrustes rotation: 0.8416
#Significance:  0.001

d.ord=capscale(mardist~1)
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

#how much variance does our model capture?
# rescaling to proportion of total variance (based on eigenvalues in ord1)
eigen.var=(ord1$CA$eig/sum(ord1$CA$eig))[names(gf$result)]
# total variance explained by model
sum(eigen.var*gf$result)

# setting the number of PCs to keep
tokeep=6
#5 with GCD, 6 with mardist

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
                    scal=12,
                    jitscale=20,
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

plot(x = xy2$lon,y=xy2$lat,col=pa2$colors,xlab=c("longitude"),ylab=c("latitude"),pch=16,xlim=c(-123.8,-122.3)) #new xlim
plot(coast$geometry, add = TRUE, col = "lightgrey")
points(x=colors2plot$lon,y=colors2plot$lat,col = "black",pch=1,cex=1.6)
points(x=colors2plot$lon,y=colors2plot$lat,col = "black",pch=1,cex=1.7)
points(x=colors2plot$lon,y=colors2plot$lat,col = "black",pch=1,cex=1.8)
points(x=colors2plot$lon,y=colors2plot$lat,col = alpha("white",0.6),pch=16,cex=1.5)
points(x=colors2plot$lon,y=colors2plot$lat,col = alpha(colors2plot$color,0.8),pch=16,cex=1.5)


