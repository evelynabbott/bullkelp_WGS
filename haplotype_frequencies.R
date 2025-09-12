

#do multiple terms at once

library(data.table)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(data.table)

#try averaging the SNPs instead
rm(list=ls())
load("GO.DF.Salinity.Rdata")
load("salinity_sig.gos.cats.Rdata")
load("lfmm64_input.Rdata")

salinity_sig.gos=salinity_sig.gos[which(salinity_sig.gos$elimKS <= 0.01),]
#salinity_sig.gos=salinity_sig.gos[salinity_sig.gos$elimKS %in% head(sort(salinity_sig.gos$elimKS), 3),]

goterm = salinity_sig.gos$GO.ID
maketab = GO.DF.SAMPLES[which(GO.DF.SAMPLES$go_acc %in% goterm),]
#write.table(maketab,file="salinity_GO2prot_supplemental.tsv",col.names = T,sep="\t")
protloc = GO.DF.SAMPLES[which(GO.DF.SAMPLES$go_acc %in% goterm),]$proteinId

snpPosFil<-hit.DF.SAMPLES[ hit.DF.SAMPLES$proteinId%in% protloc, ] #get locations of protein IDs matched with GO terms

#Scaffold and positions of the snps within the scan distance (1000bp)
SNPsterm<-unique(paste0("SCAF","_",snpPosFil$SCAFF,".",snpPosFil$SNP.POS))

#reading the raw genotype data from a different environment (created in script lfmm for Evelyn)
gen.imp<-gen.ImpN.df
#read population codes by groups
Sample.Names<-data.frame("V1"= rownames(Env.dataP))
PopCode<-substr(Sample.Names$V1,1,2)

Group1<-c("FW","MC","JT","AH","PT","KR","FB")
Group2<-c("DB","SH","PP","HI","ED","SM","EB")
Group3<-c("LP","PV","SB","DI","SQ")

GenGroupK3<-1:length(Sample.Names)
GenGroupK3[PopCode%in%Group1]<-"SJF"
GenGroupK3[PopCode%in%Group2]<-"WB"
GenGroupK3[PopCode%in%Group3]<-"SPS"

genotypesterm<-data.frame(gen.imp[,which(colnames(gen.imp)%in%SNPsterm)]) #allele per ind
#assign A or B for ref/alt
genotypesterm = data.frame(lapply(genotypesterm, function(x) {
  gsub(1,"a",x)
}))

genotypesterm = data.frame(lapply(genotypesterm, function(x) {
  gsub(0,"b",x) 
}))

rownames(genotypesterm) = rownames(gen.imp)

#get list of columns to apply haplotype paste
table(substr(colnames(genotypesterm),1,10))
repd = "SCAF_19.96"
haps = data.frame("hap" = apply( genotypesterm[ , substr(colnames(genotypesterm),1,10) %in% repd ] , 1 , paste , collapse = "" ))
unique(haps$hap)

genotypesterm <- genotypesterm %>% select(-contains(repd))
genotypesterm = cbind(genotypesterm,haps)

snpbypop=apply(genotypesterm,2,function(y){ tapply(y,GenGroupK3,function(x){data.frame(table(x)/length(x))})}) #freq of each SNP by cluster

#snpbypop[[1]][[1]]$x

for (i in seq_along(snpbypop[[1]])) {
  snpbypop[[1]][[i]]$site = sub("\\..*", "", names(snpbypop[[1]][i]))
  snpbypop[[1]][[i]]$geno = paste0(sub("\\..*", "", names(snpbypop)[1]),".",snpbypop[[1]][[i]]$x)
  snpbypop[[1]][[i]]$position = sub("\\..*", "", names(snpbypop)[1])
}

for (i in seq_along(snpbypop[[2]])) {
  snpbypop[[2]][[i]]$site = sub("\\..*", "", names(snpbypop[[2]][i]))
  snpbypop[[2]][[i]]$geno = paste0(sub("\\..*", "", names(snpbypop)[2]),".",snpbypop[[2]][[i]]$x)
  snpbypop[[2]][[i]]$position = sub("\\..*", "", names(snpbypop)[2])
}

for (i in seq_along(snpbypop[[3]])) {
  snpbypop[[3]][[i]]$site = sub("\\..*", "", names(snpbypop[[3]][i]))
  snpbypop[[3]][[i]]$geno = paste0(sub("\\..*", "", names(snpbypop)[3]),".",snpbypop[[3]][[i]]$x)
  snpbypop[[3]][[i]]$position = sub("\\..*", "", names(snpbypop)[3])
}

a=bind_rows(snpbypop[[1]])
b=bind_rows(snpbypop[[2]])
c=bind_rows(snpbypop[[3]])

pie.df=rbind(a,b,c)

pie.df$sps = ifelse(pie.df$site == "SPS" & pie.df$Freq > 0.5, "sps.max",NA)
sps.max=pie.df[which(pie.df$sps == "sps.max"),]$geno
pie.df$sps[pie.df$geno %in% sps.max] <- "sps.max"

vals = pie.df[which(pie.df$sps == "sps.max"),]$geno

scaf2 = pie.df[pie.df$position == "SCAF_2",]
test2=aggregate(Freq ~ geno, FUN = "sum", data = scaf2)
test2 = test2[order(desc(test2$Freq)),]
test2$sps = c(1:nrow(test2))
test2$sps = ifelse(test2$geno %in% vals, "sps.max",test2$sps)

scaf3 = pie.df[pie.df$position == "SCAF_3",]
test3=aggregate(Freq ~ geno, FUN = "sum", data = scaf3)
test3 = test3[order(desc(test3$Freq)),]
test3$sps = c(1:nrow(test3))
test3$sps = ifelse(test3$geno %in% vals, "sps.max",test3$sps)

scafhap = pie.df[pie.df$position == "hap",]
testhap=aggregate(Freq ~ geno, FUN = "sum", data = scafhap)
testhap = testhap[order(desc(testhap$Freq)),]
testhap$sps = c(1:nrow(testhap))
testhap$sps = ifelse(testhap$geno %in% vals, "sps.max",testhap$sps)

testall = rbind(test2,test3,testhap)
testall = testall[,c(1,3)]

pie.df = pie.df[,-c(6)]

pie.df = full_join(pie.df,testall,by="geno")

pie.df$site = factor(pie.df$site,levels = c("SPS","WB","SJF"))
pie.df$sps = factor(pie.df$sps,levels = c("sps.max",2,3,4,5))

pie.df.sal = pie.df

#save(pie.df.sal,file="pie.df.sal.Rdata")

ggplot(pie.df.sal,aes(x="",y=Freq,group=site,fill=sps))+
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y")+
  theme_classic()+
  scale_fill_manual(values=c("#FF6E6F","#00C2C3","#35B616","#FC4AB3","#CE8E1A"))+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        panel.border=element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title = element_blank(),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA),
        legend.position = "none")+
  theme(panel.spacing.x=unit(-1.5, "lines"),panel.spacing.y=unit(-.4, "lines"))+
  facet_wrap(~ position + site)

# sal.maketab = maketab
# sal.snpPosFil = snpPosFil
# save(sal.maketab,sal.snpPosFil,file="salinity_tabs4sup.Rdata")

#PAR -----
rm(list=ls())
load("lfmm64_input.Rdata")
load("GO.DF.PARsummer.Rdata")
load("PARsummer_sig.gos.cats.Rdata")

parsummer_sig.gos=parsummer_sig.gos[parsummer_sig.gos$elimKS %in% head(sort(parsummer_sig.gos$elimKS), 3),]

goterm = parsummer_sig.gos$GO.ID
maketab = GO.DF.SAMPLES[which(GO.DF.SAMPLES$go_acc %in% goterm),]
write.table(maketab,file="parsummer_GO2prot_supplemental.tsv",col.names = T,sep="\t")
protloc = GO.DF.SAMPLES[which(GO.DF.SAMPLES$go_acc %in% goterm),]$proteinId
snpPosFil<-hit.DF.SAMPLES[ hit.DF.SAMPLES$proteinId%in% protloc, ] #get locations of protein IDs matched with GO terms

#Scaffold and positions of the snps within the scan distance (1000bp)
SNPsterm<-unique(paste0("SCAF","_",snpPosFil$SCAFF,".",snpPosFil$SNP.POS))

#reading the raw genotype data from a different environment (created in script lfmm for Evelyn)
gen.imp<-gen.ImpN.df
#read population codes by groups
Sample.Names<-data.frame("V1"= rownames(Env.dataP))
PopCode<-substr(Sample.Names$V1,1,2)

Group1<-c("FW","MC","JT","AH","PT","KR","FB")
Group2<-c("DB","SH","PP","HI","ED","SM","EB")
Group3<-c("LP","PV","SB","DI","SQ")

GenGroupK3<-1:length(Sample.Names)
GenGroupK3[PopCode%in%Group1]<-"SJF"
GenGroupK3[PopCode%in%Group2]<-"WB"
GenGroupK3[PopCode%in%Group3]<-"SPS"

genotypesterm<-data.frame(gen.imp[,which(colnames(gen.imp)%in%SNPsterm)]) #allele per ind

#assign A or B for ref/alt
genotypesterm = data.frame(lapply(genotypesterm, function(x) {
  gsub(1,"a",x)
}))

genotypesterm = data.frame(lapply(genotypesterm, function(x) {
  gsub(0,"b",x) 
}))

rownames(genotypesterm) = rownames(gen.imp)


table(gsub("\\..*","",colnames(genotypesterm)))

repd = "SCAF_33"
haps33 = data.frame("hap" = apply( genotypesterm[ , gsub("\\..*","",colnames(genotypesterm)) %in% repd ] , 1 , paste , collapse = "" ))
unique(haps33$hap)
colnames(haps33) = "scaf_33.hap"

# repd = "SCAF_4"
# haps4 = data.frame("hap" = apply( genotypesterm[ , gsub("\\..*","",colnames(genotypesterm)) %in% repd ] , 1 , paste , collapse = "" ))
# unique(haps4$hap)
# colnames(haps4) = "scaf_4.hap"
# 
# repd = "SCAF_33"
# haps33 = data.frame("hap" = apply( genotypesterm[ , gsub("\\..*","",colnames(genotypesterm)) %in% repd ] , 1 , paste , collapse = "" ))
# unique(haps33$hap)
# colnames(haps33) = "scaf_33.hap"

#repd = c("SCAF_30","SCAF_4","SCAF_33")

genotypesterm <- genotypesterm %>% select(-contains(repd))
#genotypesterm = cbind(genotypesterm,haps)

genotypesterm = cbind(genotypesterm,haps33)

snpbypop=apply(genotypesterm,2,function(y){ tapply(y,GenGroupK3,function(x){data.frame(table(x)/length(x))})}) #freq of each SNP by cluster

for (i in seq_along(snpbypop[[1]])) {
  snpbypop[[1]][[i]]$site = sub("\\..*", "", names(snpbypop[[1]][i]))
  snpbypop[[1]][[i]]$geno = paste0(sub("\\..*", "", names(snpbypop)[1]),".",snpbypop[[1]][[i]]$x)
  snpbypop[[1]][[i]]$position = sub("\\..*", "", names(snpbypop)[1])
}

for (i in seq_along(snpbypop[[2]])) {
  snpbypop[[2]][[i]]$site = sub("\\..*", "", names(snpbypop[[2]][i]))
  snpbypop[[2]][[i]]$geno = paste0(sub("\\..*", "", names(snpbypop)[2]),".",snpbypop[[2]][[i]]$x)
  snpbypop[[2]][[i]]$position = sub("\\..*", "", names(snpbypop)[2])
}

for (i in seq_along(snpbypop[[3]])) {
  snpbypop[[3]][[i]]$site = sub("\\..*", "", names(snpbypop[[3]][i]))
  snpbypop[[3]][[i]]$geno = paste0(sub("\\..*", "", names(snpbypop)[3]),".",snpbypop[[3]][[i]]$x)
  snpbypop[[3]][[i]]$position = sub("\\..*", "", names(snpbypop)[3])
}

a=bind_rows(snpbypop[[1]])
b=bind_rows(snpbypop[[2]])
c=bind_rows(snpbypop[[3]])

pie.df=rbind(a,b,c)

#set factor levels for colors - greatest frequency in SPS, followed by greatest - least in others
pie.df$sps = ifelse(pie.df$site == "SPS" & pie.df$Freq > 0.5, "sps.max",NA)
sps.max=pie.df[which(pie.df$sps == "sps.max"),]$geno
pie.df$sps[pie.df$geno %in% sps.max] <- "sps.max"

vals = pie.df[which(pie.df$sps == "sps.max"),]$geno

scaf3 = pie.df[pie.df$position == "SCAF_3",]
test3=aggregate(Freq ~ geno, FUN = "sum", data = scaf3)
test3 = test3[order(desc(test3$Freq)),]
test3$sps = c(1:nrow(test3))
test3$sps = ifelse(test3$geno %in% vals, "sps.max",test3$sps)

scaf36 = pie.df[pie.df$position == "SCAF_36",]
test36=aggregate(Freq ~ geno, FUN = "sum", data = scaf36)
test36 = test36[order(desc(test36$Freq)),]
test36$sps = c(1:nrow(test36))
test36$sps = ifelse(test36$geno %in% vals, "sps.max",test36$sps)

scaf33 = pie.df[pie.df$position == "scaf_33",]
test33=aggregate(Freq ~ geno, FUN = "sum", data = scaf33)
test33 = test33[order(desc(test33$Freq)),]
test33$sps = c(1:nrow(test33))
test33$sps = ifelse(test33$geno %in% vals, "sps.max",test33$sps)

testall = rbind(test3,test36,test33)
testall = testall[,c(1,3)]

pie.df = pie.df[,-c(6)]

pie.df = full_join(pie.df,testall,by="geno")

pie.df$site = factor(pie.df$site,levels = c("SPS","WB","SJF"))
pie.df$sps = factor(pie.df$sps,levels = c("sps.max",2,3,4,5,6))

pie.df.par=pie.df
#save(pie.df.par,file="pie.df.par.Rdata")

ggplot(pie.df.par,aes(x="",y=Freq,group=site,fill=sps))+
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y")+
  theme_classic()+
  scale_fill_manual(values=c("#FF6E6F","#00C2C3","#35B616","#FC4AB3","#CE8E1A","darkblue"))+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        panel.border=element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title = element_blank(),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA),
        legend.position = "none")+
  theme(panel.spacing.x=unit(-1.5, "lines"),panel.spacing.y=unit(-.4, "lines"))+
  facet_wrap(~ position + site)

# par.maketab = maketab
# par.snpPosFil = snpPosFil
# save(par.maketab,par.snpPosFil,file="parsummer_tabs4sup.Rdata")


#kd490 -----
rm(list=ls())
load("GO.DF.kd490spring.Rdata")
#load("input_site_enrichment.Rdata")
load("kd490spring_sig.gos.cats.Rdata")
load("lfmm64_input.Rdata")

kd490spring_sig.gos=kd490spring_sig.gos[which(kd490spring_sig.gos$elimKS <= 0.01),]

goterm = kd490spring_sig.gos$GO.ID
maketab = GO.DF.SAMPLES[which(GO.DF.SAMPLES$go_acc %in% goterm),]
write.table(maketab,file="kd490spring_GO2prot_supplemental.tsv",col.names = T,sep="\t")
protloc = GO.DF.SAMPLES[which(GO.DF.SAMPLES$go_acc %in% goterm),]$proteinId
snpPosFil<-hit.DF.SAMPLES[ hit.DF.SAMPLES$proteinId%in% protloc, ] #get locations of protein IDs matched with GO terms

#Scaffold and positions of the snps within the scan distance (1000bp)
SNPsterm<-unique(paste0("SCAF","_",snpPosFil$SCAFF,".",snpPosFil$SNP.POS))

#reading the raw genotype data from a different environment (created in script lfmm for Evelyn)
gen.imp<-gen.ImpN.df
#read population codes by groups
Sample.Names<-data.frame("V1"= rownames(Env.dataP))
PopCode<-substr(Sample.Names$V1,1,2)

Group1<-c("FW","MC","JT","AH","PT","KR","FB")
Group2<-c("DB","SH","PP","HI","ED","SM","EB")
Group3<-c("LP","PV","SB","DI","SQ")

GenGroupK3<-1:length(Sample.Names)
GenGroupK3[PopCode%in%Group1]<-"SJF"
GenGroupK3[PopCode%in%Group2]<-"WB"
GenGroupK3[PopCode%in%Group3]<-"SPS"

genotypesterm<-data.frame(gen.imp[,which(colnames(gen.imp)%in%SNPsterm)]) #allele per ind

#assign A or B for ref/alt
genotypesterm = data.frame(lapply(genotypesterm, function(x) {
  gsub(1,"a",x)
}))

genotypesterm = data.frame(lapply(genotypesterm, function(x) {
  gsub(0,"b",x) 
}))

rownames(genotypesterm) = rownames(gen.imp)

snpbypop=apply(genotypesterm,2,function(y){ tapply(y,GenGroupK3,function(x){data.frame(table(x)/length(x))})}) #freq of each SNP by cluster

for (i in seq_along(snpbypop[[1]])) {
  snpbypop[[1]][[i]]$site = sub("\\..*", "", names(snpbypop[[1]][i]))
  snpbypop[[1]][[i]]$geno = paste0(sub("\\..*", "", names(snpbypop)[1]),".",snpbypop[[1]][[i]]$x)
  snpbypop[[1]][[i]]$position = sub("\\..*", "", names(snpbypop)[1])
}

for (i in seq_along(snpbypop[[2]])) {
  snpbypop[[2]][[i]]$site = sub("\\..*", "", names(snpbypop[[2]][i]))
  snpbypop[[2]][[i]]$geno = paste0(sub("\\..*", "", names(snpbypop)[2]),".",snpbypop[[2]][[i]]$x)
  snpbypop[[2]][[i]]$position = sub("\\..*", "", names(snpbypop)[2])
}

for (i in seq_along(snpbypop[[3]])) {
  snpbypop[[3]][[i]]$site = sub("\\..*", "", names(snpbypop[[3]][i]))
  snpbypop[[3]][[i]]$geno = paste0(sub("\\..*", "", names(snpbypop)[3]),".",snpbypop[[3]][[i]]$x)
  snpbypop[[3]][[i]]$position = sub("\\..*", "", names(snpbypop)[3])
}

a=bind_rows(snpbypop[[1]])
b=bind_rows(snpbypop[[2]])
c=bind_rows(snpbypop[[3]])

pie.df=rbind(a,b,c)

pie.df$sps = ifelse(pie.df$site == "SPS" & pie.df$Freq > 0.5, "sps.max",NA)
sps.max=pie.df[which(pie.df$sps == "sps.max"),]$geno
pie.df$sps[pie.df$geno %in% sps.max] <- "sps.max"

vals = pie.df[which(pie.df$sps == "sps.max"),]$geno

scaf10 = pie.df[pie.df$position == "SCAF_10",]
test10=aggregate(Freq ~ geno, FUN = "sum", data = scaf10)
test10 = test10[order(desc(test10$Freq)),]
test10$sps = c(1:nrow(test10))
test10$sps = ifelse(test10$geno %in% vals, "sps.max",test10$sps)

scaf12 = pie.df[pie.df$position == "SCAF_12",]
test12=aggregate(Freq ~ geno, FUN = "sum", data = scaf12)
test12 = test12[order(desc(test12$Freq)),]
test12$sps = c(1:nrow(test12))
test12$sps = ifelse(test12$geno %in% vals, "sps.max",test12$sps)

scaf17 = pie.df[pie.df$position == "SCAF_17",]
test17=aggregate(Freq ~ geno, FUN = "sum", data = scaf17)
test17 = test17[order(desc(test17$Freq)),]
test17$sps = c(1:nrow(test17))
test17$sps = ifelse(test17$geno %in% vals, "sps.max",test17$sps)

testall = rbind(test10,test12,test17)
testall = testall[,c(1,3)]

pie.df = pie.df[,-c(6)]

pie.df = full_join(pie.df,testall,by="geno")

pie.df$site = factor(pie.df$site,levels = c("SPS","WB","SJF"))
pie.df$sps = factor(pie.df$sps,levels = c("sps.max",2))

pie.df.kd490 = pie.df
#save(pie.df.kd490,file="pie.df.kd490.Rdata")

ggplot(pie.df.kd490,aes(x="",y=Freq,group=site,fill=sps))+
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y")+
  theme_classic()+
  scale_fill_manual(values=c("#FF6E6F","#00C2C3"))+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        panel.border=element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title = element_blank(),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA),
        legend.position = "none")+
  theme(panel.spacing.x=unit(-1.5, "lines"),panel.spacing.y=unit(-.4, "lines"))+
  facet_wrap(~ position + site)

# kds.maketab = maketab
# kds.snpPosFil = snpPosFil
# save(kds.maketab,kds.snpPosFil,file="kd490spring_tabs4sup.Rdata")



#plot together ----------
rm(list=ls())

load("pie.df.sal.Rdata")
load("pie.df.par.Rdata")
load("pie.df.kd490.Rdata")


sal.plot=ggplot(pie.df.sal,aes(x="",y=Freq,group=site,fill=sps))+
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y")+
  theme_classic()+
  scale_fill_manual(values=c("#FF6E6F","#00C2C3","#35B616","#FC4AB3","#CE8E1A"))+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        panel.border=element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title = element_blank(),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA),
        legend.position = "none")+
  theme(panel.spacing.x=unit(-1.5, "lines"),panel.spacing.y=unit(-.4, "lines"))+
  facet_wrap(~ position + site)


par.plot=ggplot(pie.df.par,aes(x="",y=Freq,group=site,fill=sps))+
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y")+
  theme_classic()+
  scale_fill_manual(values=c("#FF6E6F","#00C2C3","#35B616","#FC4AB3","#CE8E1A","darkblue","#8F6FFD","yellow","#96A418","orchid4"))+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        panel.border=element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title = element_blank(),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA),
        legend.position = "none")+
  theme(panel.spacing.x=unit(-1.5, "lines"),panel.spacing.y=unit(-.4, "lines"))+
  facet_wrap(~ position + site)

kd490.plot=ggplot(pie.df.kd490,aes(x="",y=Freq,group=site,fill=sps))+
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y")+
  theme_classic()+
  scale_fill_manual(values=c("#FF6E6F","#00C2C3"))+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        panel.border=element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title = element_blank(),
        panel.background = element_rect(fill='transparent'), #transparent panel bg
        plot.background = element_rect(fill='transparent', color=NA),
        legend.position = "none")+
  theme(panel.spacing.x=unit(-1.5, "lines"),panel.spacing.y=unit(-.4, "lines"))+
  facet_wrap(~ position + site)

grid.arrange(sal.plot,par.plot,kd490.plot,nrow=3)






