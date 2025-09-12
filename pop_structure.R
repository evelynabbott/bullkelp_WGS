#PCA and admix

library(DESeq2)
library(tidyverse)
library(ggplot2)
library(vegan)
library(scatterpie)

rm(list=ls())
#start with PCA

load("samples65.Rdata")

#load the data - ibs matrix and list of bams ---------------------------
bams_tacc=read.table("bams.nr")
#rename the replicate
bams_tacc$V1=str_replace(bams_tacc$V1,"JT08FA6_S42_sort1.rmd.bam","JT08FA6.rep_S42_sort1.rmd.bam")
bams_tacc$V1=gsub('_S.*','',bams_tacc$V1) #rename to match metadata
length(unique(bams_tacc$V1)) #make sure all names are unique

ibsMat=read.table("myresult.ibsMat")
colnames(ibsMat) = bams_tacc$V1
rownames(ibsMat) = bams_tacc$V1

load("metadata.Rdata")
meta$ind = colnames(ibsMat)

#remove bad samples
ibsMat = ibsMat[which(rownames(ibsMat) %in% samples65),which(colnames(ibsMat) %in% samples65)]
meta = meta[which(meta$ind %in% samples65),]

#add clusters to metadata
Group1<-c("FW","MC","JT","AH","PT","KR","FB")
Group2<-c("DB","SH","PP","HI","ED","SM","EB")
Group3<-c("LP","PV","SB","DI","SQ")

which(meta$site%in%Group1)
which(meta$site%in%Group2)
which(meta$site%in%Group3)

meta$clust = ifelse(meta$site%in%Group1,1,
                    ifelse(meta$site%in%Group2,2,3))


#PCA run capscale --------
pp0_ad=capscale(ibsMat~1)
plot(pp0_ad$CA$eig/sum(pp0_ad$CA$eig))
scores=pp0_ad$CA$u

#lab=rownames(ibsMat)
lab=meta$clust

#inital plot to identify outliers
axes2plot=c(1,2)
plot(scores[,axes2plot],pch=19,main="",sub="",cex=0.2)
#ordispider(scores[,axes2plot],lab,label=T,cex=0.5)
ordihull(scores[,axes2plot],lab,label=T,cex=0.8)

pcs=as.data.frame(scores)
pcs = cbind(meta,pcs)
pcs$clust=as.factor(pcs$clust)

#ggplot
ggplot(pcs,aes(x=MDS1,y=MDS2,color=clust))+
  geom_point(size=5,alpha=0.7)+
  theme_classic()+
  #coord_fixed()+
  theme(legend.position = "none")+
  scale_color_manual(values=c("hotpink","orchid4","darkorange3"))


#admixture ----------------
# rm(list=ls())
# source("plot_admixture_v4_function.R")
# load("metadata.Rdata")
# npops <- 3 #choose number of populations
# 
# inName <- paste0('pops_k', npops, '.qopt') #inName corresponds to .qopt files to upload later
# #pops
# inds=read.table("bams.nr",sep="\t") #bams used
# inds=data.frame(lapply(inds, gsub, pattern="_S.*",replacement = "")) #make file names match metadata
# 
# pops=data.frame("id"=meta$ind,"habitat"=meta$region)
# write.table(pops, file = "pops",sep = "\t",col.names = c("id","habitat")) #save as tab delimited
# pops=read.table("pops") #read as tab delimited table
# 
# 
# npops=as.numeric(sub("\\D+(\\d+)\\..+","\\1",inName))
# tbl=read.table(paste(inName,sep=""),header=F)
# 
# i2p=read.table("pops")
# names(i2p)=c("ind","pop")
# tbl=cbind(tbl,i2p)
# rownames(tbl) <- tbl$ind <- sub("(.*?)\\..*$", "\\1", tbl$ind)
# head(tbl,20) # this is how the resulting dataset must look
# 
# #save for pies
# sites=meta$site
# pies = cbind(tbl,sites)
# pies = pies[-c(1),]
# save(pies, file = "k3_pies_for_maps.Rdata")
# 
# #plot differently
# tbl$pop=factor(tbl$pop,levels=c("PS","U","SJF"))
# #tbl$pop=factor(tbl$pop)
# 
# plotAdmixture(data=tbl,npops=npops,grouping.method="distance") +
#   scale_fill_manual(values = c("hotpink2","darkorange3","orchid4"))


#make pies --------------------------------
rm(list=ls())
load("k3_pies_for_maps.Rdata")
load("samples65.Rdata")

pies = pies[which(rownames(pies) %in% samples65),]

pie_df = pies

#AH
pie = subset(pie_df,pie_df$site=="AH")
pie = pie[,c(1:3)]
AH = data.frame(site="AH",t(colSums(pie)))

#DB
pie = subset(pie_df,pie_df$site=="DB")
pie = pie[,c(1:3)]
DB = data.frame(site="DB",t(colSums(pie)))

#DI
pie = subset(pie_df,pie_df$site=="DI")
pie = pie[,c(1:3)]
DI = data.frame(site="DI",t(colSums(pie)))

#EB
pie = subset(pie_df,pie_df$site=="EB")
pie = pie[,c(1:3)]
EB = data.frame(site="EB",t(colSums(pie)))

#ED
pie = subset(pie_df,pie_df$site=="ED")
pie = pie[,c(1:3)]
ED = data.frame(site="ED",t(colSums(pie)))

#FB
pie = subset(pie_df,pie_df$site=="FB")
pie = pie[,c(1:3)]
FB = data.frame(site="FB",t(colSums(pie)))

#FW
pie = subset(pie_df,pie_df$site=="FW")
pie = pie[,c(1:3)]
FW = data.frame(site="FW",t(colSums(pie)))

#HI
pie = subset(pie_df,pie_df$site=="HI")
pie = pie[,c(1:3)]
HI = data.frame(site="HI",t(colSums(pie)))

#JT
pie = subset(pie_df,pie_df$site=="JT")
pie = pie[,c(1:3)]
JT = data.frame(site="JT",t(colSums(pie)))

#KR
pie = subset(pie_df,pie_df$site=="KR")
pie = pie[,c(1:3)]
KR = data.frame(site="KR",t(colSums(pie)))

#LP
pie = subset(pie_df,pie_df$site=="LP")
pie = pie[,c(1:3)]
LP = data.frame(site="LP",t(colSums(pie)))

#MC
pie = subset(pie_df,pie_df$site=="MC")
pie = pie[,c(1:3)]
MC = data.frame(site="MC",t(colSums(pie)))

#PP
pie = subset(pie_df,pie_df$site=="PP")
pie = pie[,c(1:3)]
PP = data.frame(site="PP",t(colSums(pie)))

#PT
pie = subset(pie_df,pie_df$site=="PT")
pie = pie[,c(1:3)]
PT = data.frame(site="PT",t(colSums(pie)))

#PV
pie = subset(pie_df,pie_df$site=="PV")
pie = pie[,c(1:3)]
PV = data.frame(site="PV",t(colSums(pie)))

#SB
pie = subset(pie_df,pie_df$site=="SB")
pie = pie[,c(1:3)]
SB = data.frame(site="SB",t(colSums(pie)))

#SH
pie = subset(pie_df,pie_df$site=="SH")
pie = pie[,c(1:3)]
SH = data.frame(site="SH",t(colSums(pie)))

#SM
pie = subset(pie_df,pie_df$site=="SM")
pie = pie[,c(1:3)]
SM = data.frame(site="SM",t(colSums(pie)))

#SQ
pie = subset(pie_df,pie_df$site=="SQ")
pie = pie[,c(1:3)]
SQX = data.frame(site="SQX",t(colSums(pie)))


site_pies=rbind(AH,DB,DI,EB,ED,FB,FW,HI,JT,KR,LP,MC,PP,PT,PV,SB,SH,SM,SQX)

sites<-read.delim("Sample coordinates.txt")
colnames(sites) = c("name","site","lat","long")

keeps = site_pies$site
sites1 = sites[sites$site %in% keeps,]

setdiff(site_pies$site,sites1$site) #0

coords = sites1 %>% select(site,lat,long)

pie_gg=full_join(coords,site_pies,by="site")
pie_gg$radius = rep(1,19)

test_pies = pie_gg %>% select(2:6)
colnames(test_pies) = c("y","x","a","b","c")


p=ggplot()+
  geom_scatterpie(aes(x=x,y=y),
                  data = test_pies,
                  cols=c("a","b","c"),
                  pie_scale = 1.2,
                  alpha = 0.8)+
  theme(legend.position="none")+
  scale_fill_manual(values = c("orchid4","darkorange3","hotpink2"))+
  coord_fixed()+
  theme(
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent') #transparent legend panel
  )


ggsave(
  plot = p,
  filename = "tr_tst3.png",
  bg = "transparent"
)


#heterozygosity
het = read.table("het_vals.txt",header = T)
colnames(het) = c("site","het")


stats = full_join(pie_gg,het, by="site")

ggplot(stats,aes(long,lat,color=het))+
  geom_point(size=5)+
  scale_color_viridis()





