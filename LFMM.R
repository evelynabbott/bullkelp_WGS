


#LFMM 65 samples


#plotting RDA, LFMM, and identifying SNPs
library(devtools)
library(ggvegan)
library(lfmm)
library(qvalue)
library(stringr)
library(dplyr)
library(vegan)
library(ggplot2)
library(qqman)

rm(list=ls())

#use samples with >=7x coverage
#load("rda_lfmm_files.Rdata")
load("lfmm64_input.Rdata")
#load("env65.noRast.Rdata")

pairs(Env.dataP)

#keep important
Env.dataP = Env.dataP[,-c(4,7,8,9)]


#for REHH, get rid of sps
# SPS<-c("LP","PV","SB","DI","SQ")
# Env.dataP$pop = substr(rownames(Env.dataP),1,2)
# gen.ImpN.df$pop = substr(rownames(gen.ImpN.df),1,2)
# `%nin%` <- Negate(`%in%`)
# Env.dataP = Env.dataP[Env.dataP$pop %nin% SPS,]
# Env.dataP = Env.dataP[,-c(8)]
# gen.ImpN.df = gen.ImpN.df[gen.ImpN.df$pop %nin% SPS,]
# gen.ImpN.df = gen.ImpN.df[,-c(712251)]


pred.pcaNereo <- rda(Env.dataP, scale=T)
summary(pred.pcaNereo)$cont
screeplot(pred.pcaNereo,  bstick=TRUE, main = "Screeplot: Eigenvalues of Puget Soun Env Variables 3")

# 1 2 3
## correlations between the PC axis and predictors:----
Corr.PCs.and.EnvVar<-round(scores(pred.pcaNereo, choices=1:3, display="species", scaling=0), digits=3)

#                   PC1    PC2    PC3
# kd490.Winter      0.516 -0.146  0.203
# PAR.Summer2       0.116 -0.422 -0.308
# kd490.Spring      0.391  0.282  0.312
# PAR.Winter2      -0.260 -0.498 -0.003
# SalinityAVG2     -0.306 -0.253  0.494
# SST.month.Winter  0.317 -0.459 -0.120
# SST.August        0.529 -0.094 -0.189
# ChlA.Winter       0.159 -0.063  0.650
# ChlA.Spring       0.027  0.431 -0.226

K=3

screeplot(gen.pca, main = "Screeplot of Genetic Data with Broken Stick", bstick=TRUE, type="barplot")

#We’ll store our synthetic PC axis predictor as pred.PC1 for use in LFMM.
pred.PC1 <- scores(pred.pcaNereo, choices=1, display="sites", scaling=0)
pred.PC2 <- scores(pred.pcaNereo, choices=2, display="sites", scaling=0)
#pred.PC3 <- scores(pred.pcaNereo, choices=3, display="sites", scaling=0)
# pred.PC4 <- scores(pred.pcaNereo, choices=4, display="sites", scaling=0)

#pred.PC1 = as.numeric(pred.PC1)

Nereo.lfmmPC1.k3 <- lfmm_ridge(Y=gen.ImpN.df, X=pred.PC1, K=3) ## c.ange K as you see fit
Nereo.lfmmPC2.k3 <- lfmm_ridge(Y=gen.ImpN.df, X=pred.PC2, K=3) ## c.ange K as you see fit
#Nereo.lfmmPC3.k3 <- lfmm_ridge(Y=gen.ImpN.df, X=pred.PC3, K=3)
# Nereo.lfmmPC3.k4 <- lfmm_ridge(Y=gen.ImpN.df, X=pred.PC4, K=3)

#LFMM per variable
# Nereo.lfmmSST.August.k3 <- lfmm_ridge(Y=gen.ImpN.df, X=Env.dataP$SST.August, K=3) ## c.ange K as you see fit
Nereo.lfmmSST.Salinity.k3 <- lfmm_ridge(Y=gen.ImpN.df, X=Env.dataP$SalinityAVG2, K=3)
Nereo.lfmmSST.Winter.k3 <- lfmm_ridge(Y=gen.ImpN.df, X=Env.dataP$SST.month.Winter, K=3)
Nereo.lfmm.Parsummer.k3 <- lfmm_ridge(Y=gen.ImpN.df, X=Env.dataP$PAR.Summer2, K=3)
# Nereo.lfmm.ChlAspring.k3 <- lfmm_ridge(Y=gen.ImpN.df, X=Env.dataP$ChlA.Spring, K=3)
# Nereo.lfmm.ChlAwinter.k3 <- lfmm_ridge(Y=gen.ImpN.df, X=Env.dataP$ChlA.Winter, K=3)
Nereo.lfmm.kd490winter.k3 <- lfmm_ridge(Y=gen.ImpN.df, X=Env.dataP$kd490.Winter, K=3)
Nereo.lfmm.kd490spring.k3 <- lfmm_ridge(Y=gen.ImpN.df, X=Env.dataP$kd490.Spring, K=3)
# Nereo.lfmm.Parwinter.k3 <- lfmm_ridge(Y=gen.ImpN.df, X=Env.dataP$PAR.Winter2, K=3)

#PCs pvalues
Nereo.pv.pc1.k3 <- lfmm_test(Y=gen.ImpN.df, X=pred.PC1, lfmm=Nereo.lfmmPC1.k3, calibrate="gif")
Nereo.pv.pc2.k3 <- lfmm_test(Y=gen.ImpN.df, X=pred.PC2, lfmm=Nereo.lfmmPC2.k3, calibrate="gif")
#Nereo.pv.pc3.k3<- lfmm_test(Y=gen.ImpN.df, X=pred.PC3, lfmm=Nereo.lfmmPC3.k3, calibrate="gif")

#variable pvals
Nereo.pv.Salinity.k3<- lfmm_test(Y=gen.ImpN.df, X=Env.dataP$SalinityAVG, lfmm=Nereo.lfmmSST.Salinity.k3, calibrate="gif")
#Nereo.pv.AugustSST.k3<- lfmm_test(Y=gen.ImpN.df, X=Env.dataP$SST.August, lfmm=Nereo.lfmmSST.August.k3, calibrate="gif")
Nereo.pv.WinterSST.k3<- lfmm_test(Y=gen.ImpN.df, X=Env.dataP$SST.month.Winter, lfmm=Nereo.lfmmSST.Winter.k3, calibrate="median+MAD")
Nereo.pv.Parsummer.k3<- lfmm_test(Y=gen.ImpN.df, X=Env.dataP$PAR.Summer2, lfmm=Nereo.lfmm.Parsummer.k3, calibrate="gif")
#Nereo.pv.ChlAspring.k3<- lfmm_test(Y=gen.ImpN.df, X=Env.dataP$ChlA.Spring, lfmm=Nereo.lfmm.ChlAspring.k3, calibrate="gif")
#Nereo.pv.ChlAwinter.k3<- lfmm_test(Y=gen.ImpN.df, X=Env.dataP$ChlA.Winter, lfmm=Nereo.lfmm.ChlAwinter.k3, calibrate="gif")
Nereo.pv.kd490winter.k3<- lfmm_test(Y=gen.ImpN.df, X=Env.dataP$kd490.Winter, lfmm=Nereo.lfmm.kd490winter.k3, calibrate="gif")
Nereo.pv.kd490spring.k3<- lfmm_test(Y=gen.ImpN.df, X=Env.dataP$kd490.Spring, lfmm=Nereo.lfmm.kd490spring.k3, calibrate="gif")
# Nereo.pv.Parwinter.k3<- lfmm_test(Y=gen.ImpN.df, X=Env.dataP$PAR.Winter2, lfmm=Nereo.lfmm.Parwinter.k3, calibrate="gif")



# #names(Nereo.pv) # this object includes raw z-scores and p-values, as well as GIF-calibrated scores and p-values
# #GIF values----
Nereo.pv.pc1.k3$gif
# #1.334212 
Nereo.pv.pc2.k3$gif
# #2.076674

Nereo.pv.Salinity.k3$gif
#[1] 2.514396
Nereo.pv.WinterSST.k3$gif
#[1] 0.4332478 **very low
Nereo.pv.Parsummer.k3$gif
#[1] 1.865715
Nereo.pv.kd490winter.k3$gif
#[1] 1.304167
Nereo.pv.kd490spring.k3$gif
#1.272476

#An appropriately calibrated set of tests will have a GIF of around 1. 
#The elevated GIF for our tests indicates that the results may be overly liberal
#in identifying candidate SNPs. If the GIF is less than one, the test may be too conservative.

#For .K3
# x11()
# par(mfrow=c(2,1))
# hist(Nereo.pv.pc1.k3$pvalue[,1], main="Unadjusted p-values pc1 k3")
# hist(Nereo.pv.pc1.k3$calibrated.pvalue[,1], main="GIF-adjusted p-values pc1 k3")  #Some peak here, also the lowest GIF at 1.2
# dev.off()
# 
# par(mfrow=c(2,1))
# hist(Nereo.pv.pc2.k3$pvalue[,1], main="Unadjusted p-values pc2 k3")
# hist(Nereo.pv.pc2.k3$calibrated.pvalue[,1], main="GIF-adjusted p-values pc2 k3")
# 
# # x11()
# # par(mfrow=c(2,1))
# hist(Nereo.pv.pc3.k3$pvalue[,1], main="Unadjusted p-values pc3 k3")
# hist(Nereo.pv.pc3.k3$calibrated.pvalue[,1], main="GIF-adjusted p-values pc3 k3")

#For Salinity and SST August separately
# hist(Nereo.pv.AugustSST.k3$pvalue[,1], main="SST AUG Unadjusted p-values pc4 k3")
# hist(Nereo.pv.AugustSST.k3$calibrated.pvalue[,1], main="SST AUG GIF-adjusted p-values pc4 k3")
# 
# #
# hist(Nereo.pv.WinterSST.k3$pvalue[,1], main="SST WINTER Unadjusted p-values pc4 k3")
# hist(Nereo.pv.WinterSST.k3$calibrated.pvalue[,1], main="SST WINTER GIF-adjusted p-values pc4 k3")
# 
# hist(Nereo.pv.Salinity.k3$pvalue[,1], main="Salinity Unadjusted p-values pc4 k3")
# hist(Nereo.pv.Salinity.k3$calibrated.pvalue[,1], main="Salinity GIF-adjusted p-values pc4 k3")
# 
# hist(Nereo.pv.Salinity.k3$pvalue[,1], main="Salinity Unadjusted p-values pc4 k3")
# hist(Nereo.pv.Salinity.k3$calibrated.pvalue[,1], main="Salinity GIF-adjusted p-values pc4 k3")



# 
#q vals ------
#K3
Nereo.qv.pc1.k3 <- qvalue(Nereo.pv.pc1.k3$calibrated.pvalue)$qvalues
Nereo.qv.pc2.k3 <- qvalue(Nereo.pv.pc2.k3$calibrated.pvalue)$qvalues
#Nereo.qv.pc3.k3 <- qvalue(Nereo.pv.pc3.k3$calibrated.pvalue)$qvalues

# #Nereo.qv.sst.August.k3 <- qvalue(Nereo.pv.AugustSST.k3$calibrated.pvalue)$qvalues
Nereo.qv.salinity.k3 <- qvalue(Nereo.pv.Salinity.k3$calibrated.pvalue)$qvalues
Nereo.qv.WinterSST.k3 <- qvalue(Nereo.pv.WinterSST.k3$calibrated.pvalue)$qvalues
Nereo.qv.Parsummer.k3 <- qvalue(Nereo.pv.Parsummer.k3$calibrated.pvalue)$qvalues
#Nereo.qv.ChlAspring.k3 <- qvalue(Nereo.pv.ChlAspring.k3$calibrated.pvalue)$qvalues
#Nereo.qv.ChlAwinter.k3 <- qvalue(Nereo.pv.ChlAwinter.k3$calibrated.pvalue)$qvalues
Nereo.qv.kd490winter.k3 <- qvalue(Nereo.pv.kd490winter.k3$calibrated.pvalue)$qvalues
Nereo.qv.kd490spring.k3 <- qvalue(Nereo.pv.kd490spring.k3$calibrated.pvalue)$qvalues
# Nereo.qv.Parwinter.k3 <- qvalue(Nereo.pv.Parwinter.k3$calibrated.pvalue)$qvalues

#save these files to clear up space in env
save(list = ls(pattern = '.pv.'), file = "pvals4go.RData")
save(list = ls(pattern = '.qv.'), file = "qvals4go.RData")


#These can be used to extract the loci under section associated with different 
#PC1
## how many SNPs have an FDR < 10%?
length(which(Nereo.qv.pc1.k3 < 0.1)) # 1065
length(which(Nereo.qv.pc2.k3 < 0.1)) # 221
# length(which(Nereo.qv.pc3.k3 < 0.1)) # 17

length(which(Nereo.qv.WinterSST.k3 < 0.05)) # 12249
length(which(Nereo.qv.salinity.k3 < 0.05)) # 141 *** 1
length(which(Nereo.qv.kd490winter.k3 < 0.05)) # 413
length(which(Nereo.qv.kd490spring.k3 < 0.05)) # 48 * 3
length(which(Nereo.qv.Parsummer.k3 < 0.05)) # 1725 ** 2

#get stats ----------------------------------------------------------
rm(list=ls())
load("pvals4go.RData")
load("qvals4go.RData")

#scaffolds?
qv.list = mget(ls(pattern = '.qv.'))
qv.list = qv.list[c(1,2,3,6,7)]
qv.dfs = list()

for (i in 1:length(qv.list)) {
  qv.dfs[[i]] = data.frame("scaff" = gsub("\\..*","",rownames(qv.list[[i]])),
                           "pos" = gsub("^.*\\.","",rownames(qv.list[[i]])),
                           "qval" = qv.list[[i]][1:nrow(qv.list[[i]])],
                           "var" = gsub("Nereo.qv.","",names(qv.list)[i]))
  qv.dfs[[i]] = qv.dfs[[i]][which(qv.dfs[[i]][,3] <= 0.05),]
}

table(qv.dfs[[1]][,1]) #scaf 25, 19x
table(qv.dfs[[2]][,1]) #scaf 1, 57x; scaf 5, 55x; scaf 7, 47x
max(table(qv.dfs[[3]][,1])) #scaf 10 (187); 13(125); 19(121); 2(120); 22(125); 30(147); 33(170); 4(198); 5(104); 6(108);8(100)
max(table(qv.dfs[[4]][,1])) #even
max(table(qv.dfs[[5]][,1])) #even

length(unique(qv.dfs[[1]][,1])) #24 scaffolds, most on 19
length(unique(qv.dfs[[2]][,1])) #39 scaffolds, most on 1
length(unique(qv.dfs[[3]][,1])) #50 scaffolds, most on 4
length(unique(qv.dfs[[5]][,1])) #32 scaffolds, most on 4


all.q = bind_rows(qv.dfs)
all.q$scaf.num = as.numeric(gsub("SCAF_","",all.q$scaff))
all.trim = all.q[!duplicated(all.q$qval),]
all.dups = all.q[duplicated(all.q$qval),]

table(all.dups$var)

#duplicated pvals?
pv.list = mget(ls(pattern = '.pv.'))
pv.dfs <- list()
for (i in 1:length(pv.list)) {
  pv.dfs[[i]] = data.frame("snp" = rownames(pv.list[[i]]$B),
                           "pval" = pv.list[[i]]$calibrated.pvalue,
                           "var" = names(pv.list[i]))
  pv.dfs[[i]] = pv.dfs[[i]][which(pv.dfs[[i]][,2] <= 0.05),]
}

cnms = c("snp","pval","var")
pv.dfs=lapply(pv.dfs, setNames, cnms)

all.p = bind_rows(pv.dfs)
p.dups = all.p[duplicated(all.p$pval),]
length(which(paste0(all.dups$scaff,".",all.dups$pos) %in% p.dups$snp)) #6846 of 8392

table(all.p$var)

#calibrate another way?











tail(names(sort(table(all.q$scaff))), 1)
a=as.data.frame(table(all.q$scaff))
a$Freq = as.numeric(a$Freq)
big = a[a$Freq > 226,]
sum(big$Freq)

#average pvals 

#manhattan combined
qv.list = mget(ls(pattern = '.qv.'))
qv.list = qv.list[c(1,2,3,6,7)]
qv.mann = list()

for (i in 1:length(qv.list)) {
  qv.mann[[i]] = data.frame("pos" = gsub("SCAF_","",rownames(qv.list[[i]])),
                           "scaf" = gsub("^.*\\.","",rownames(qv.list[[i]])),
                           "qval" = qv.list[[i]][1:nrow(qv.list[[i]])],
                           "var" = gsub("Nereo.qv.","",names(qv.list)[i]))
}

df4mann = bind_rows(qv.mann)

df4mann$scaf = as.numeric(gsub("\\..*","",df4mann$pos))

df4mann$color = ifelse(df4mann$qval <= 0.01,"red",
                       ifelse(df4mann$scaf %% 2 == 1, "#3C3C3C",
                              ifelse(df4mann$scaf %% 2 == 0, "#ABABAB",NA)))

plot(df4mann$pos,-log(df4mann$qval))

#manhat
mann = data.frame("CHR" = df4mann$scaf,"BP" = as.numeric(gsub("^.*\\.","",df4mann$pos)),"P" = df4mann$qval)
mann$SNP = paste0("SCAF",mann$CHR,"_",mann$BP)

snpsinterest = mann[(which(mann$P < 0.05)),]
snpsinterest = snpsinterest$SNP

manhattan(mann,
          cex=1,
          suggestiveline = F,
          genomewideline = F,
          highlight = snpsinterest)

#remove snps with identical pvals


#The I copied the file sizes.genome to this working directory
sizes.genome<-read.table("sizes.genome",header=F)

# #PAR --------------
# pvalues<-Nereo.pv.Parsummer.k3$calibrated.pvalue
# qvalues<-Nereo.qv.Parsummer.k3
# 
# SCAFF<-as.numeric(str_extract(rownames(pvalues)[-1],regex("(?<=(_))[0-9]+")))
# POS <-as.numeric(str_extract(rownames(pvalues)[-1],regex("(?<=(-))[0-9]+")))
# pvalues.DF<-data.frame(SCAFF,POS,P.VAL=as.numeric(pvalues)[-1],Q.VAL=qvalues[-1])
# pvalues.DF<-data.frame(SCAFF,POS,P.VAL=as.numeric(pvalues)[-1],Q.VAL=qvalues[-1])
# SCAFF.GEN<-as.numeric(str_extract(sizes.genome[,1],regex("(?<=(_))[0-9]+")))
# sizes.genome.DF<-data.frame(SCAFF.GEN,SIZE=sizes.genome[,2])
# U.scaff.codesPvalue<-unique(pvalues.DF$SCAFF)
# n.scaffs.in.pval<-length(U.scaff.codesPvalue)
# all(U.scaff.codesPvalue%in%sizes.genome.DF$SCAFF.GEN)
# 
# pvals_allgenome.par = pvalues.DF
# pvals_allgenome.par$sig = ifelse(pvals_allgenome.par$Q.VAL > 0.1,"not_sig","sig")
# pvals_allgenome.par$color = ifelse(pvals_allgenome.par$Q.VAL > 0.1,"grey96","red")
# pvals_allgenome.par$order = c(1:1896676)
# 
# pvals_allgenome.par$scaffcol = ifelse(pvals_allgenome.par$sig == "sig", "sig",
#                                       ifelse(pvals_allgenome.par$SCAFF %% 2 == 0, "even","odd"))
# 
# 
# a=ggplot(pvals_allgenome.par,aes(order,abs(log(Q.VAL)),color=scaffcol))+
#   geom_point(alpha=0.5)+
#   scale_color_manual(values=c("grey","darkgrey","red"))+
#   theme_classic()+
#   theme(legend.position="none")
# 
# 
# #kd490 -------------
# pvalues<-Nereo.pv.kd490winter.k3$calibrated.pvalue
# qvalues<-Nereo.qv.kd490winter.k3
# 
# SCAFF<-as.numeric(str_extract(rownames(pvalues)[-1],regex("(?<=(_))[0-9]+")))
# POS <-as.numeric(str_extract(rownames(pvalues)[-1],regex("(?<=(-))[0-9]+")))
# pvalues.DF<-data.frame(SCAFF,POS,P.VAL=as.numeric(pvalues)[-1],Q.VAL=qvalues[-1])
# pvalues.DF<-data.frame(SCAFF,POS,P.VAL=as.numeric(pvalues)[-1],Q.VAL=qvalues[-1])
# SCAFF.GEN<-as.numeric(str_extract(sizes.genome[,1],regex("(?<=(_))[0-9]+")))
# sizes.genome.DF<-data.frame(SCAFF.GEN,SIZE=sizes.genome[,2])
# U.scaff.codesPvalue<-unique(pvalues.DF$SCAFF)
# n.scaffs.in.pval<-length(U.scaff.codesPvalue)
# all(U.scaff.codesPvalue%in%sizes.genome.DF$SCAFF.GEN)
# 
# pvals_allgenome.kd = pvalues.DF
# pvals_allgenome.kd$sig = ifelse(pvals_allgenome.kd$Q.VAL > 0.1,"not_sig","sig")
# pvals_allgenome.kd$color = ifelse(pvals_allgenome.kd$Q.VAL > 0.1,"grey96","red")
# pvals_allgenome.kd$order = c(1:1896676)
# 
# pvals_allgenome.kd$scaffcol = ifelse(pvals_allgenome.kd$sig == "sig", "sig",
#                                      ifelse(pvals_allgenome.kd$SCAFF %% 2 == 0, "even","odd"))
# 
# 
# ggplot(pvals_allgenome.kd,aes(order,abs(log(Q.VAL)),color=scaffcol))+
#   geom_point(alpha=0.5)+
#   scale_color_manual(values=c("grey","darkgrey","red"))+
#   theme_classic()+
#   #scale_y_continuous(limits = c(0, 50))+
#   theme(legend.position="none")
# 
# grid.arrange(a,b,nrow=2)
# 
# 
# #salinity -------------------
# pvals_allgenome.sal = pvalues.DF
# pvals_allgenome.sal$sig = ifelse(pvals_allgenome.sal$Q.VAL > 0.1,"not_sig","sig")
# pvals_allgenome.sal$color = ifelse(pvals_allgenome.sal$Q.VAL > 0.1,"grey96","red")
# pvals_allgenome.sal$order = c(1:1896676)
# pvals_allgenome.sal$scaffcol = ifelse(pvals_allgenome.sal$sig == "sig", "sig",
#                                       ifelse(pvals_allgenome.sal$SCAFF %% 2 == 0, "even","odd"))
# 
# 
# ggplot(pvals_allgenome.sal,aes(order,abs(log(Q.VAL)),color=scaffcol))+
#   geom_point(alpha=0.5)+
#   scale_color_manual(values=c("grey","darkgrey","red"))+
#   theme_classic()+
#   theme(legend.position="none")
# 
# #grid.arrange(a,b,c,d,e,nrow=3)
# 
# # df = data.frame("x" = c(1:3), "y" = c(4,5,6))
# # ggplot(pv,aes(x,y))+
# #   geom_point()
# 









# #RDA--------------------------------
# #use samples with >=7x coverage
# # load("rda_lfmm_files.Rdata")
# # save(Env.data,gen.ImpN.df,file="rda_input.Rdata")
# rm(list=ls())
# load("rda.input77.Rdata")
# 
# rdanames = colnames(gen.ImpN.df)
# save(rdanames,file="rdanames.Rdata")
# 
# #Env.dataP = Env.dataP[,c(1,4:8,10,11)] #remove high vif variables
# #remove only winter SST
# Env.dataP = Env.dataP[,c(1:8)]
# 
# nereo.rda = rda(gen.ImpN.df ~ ., data=Env.dataP, scale=T)
# nereo.rda
# save(nereo.rda,file="nereo.rda77.Rdata")
# 
# rm(list=ls())
# load("nereo.rda77.Rdata")
# load("rdanames.Rdata")
# 
# #adjust r2
# RsquareAdj(nereo.rda)
# #[1] 0.18821 - not bad including all
# #[1] 0.200316 - with high VIF removed
# 
# summary(nereo.rda)$concont
# screeplot(nereo.rda) #not much past 2
# #anova.cca(nereo.rda, by="axis")
# 
# #check variance inflation factors. Ideally most are between 5 and 10
# barplot(vif.cca(nereo.rda))
# # SalinityAVG2   SST.month.Winter       SST.August      PAR.Summer2      ChlA.Winter      ChlA.Spring 
# # 4.496628        61.137395.            52.578821       2.949314         2.677271         2.116838 
# 
# # kd490.Winter     kd490.Spring         temp_range      lat              lon 
# # 9.147727         1.607788             NA              21.441533        24.028160 
# 
# # #without sst ****
# # SalinityAVG2  PAR.Summer2  ChlA.Winter  ChlA.Spring  kd490.Winter   kd490.Spring      lat          lon  ***looks good
# # 3.393005      2.331678     1.332734     1.738064     2.198203       1.588823          3.210189     3.129768
# 
# #make plot for methods
# vif.df = data.frame("vif" = vif.cca(nereo.rda))
# vif.df$var = rownames(vif.df)
# 
# vif.df=vif.df[-c(9),]
# 
# ggplot(vif.df,aes(x=var,y=vif))+
#   geom_bar(stat="identity")
# 
# vif.df$group = ifelse(grepl("SST",vif.df$var),'SST',
#                       ifelse(grepl("kd490",vif.df$var),'kd490',
#                              ifelse(grepl("PAR",vif.df$var),'PAR',
#                                     ifelse(grepl("Chl",vif.df$var),'ChlA',
#                                            ifelse(grepl("Sal",vif.df$var),'salinity',"coordinate")))))
# 
# #vif.df$group = factor(vif.df$group,levels=c("SST","latlon","kd490","salinity","PAR","ChlA"))
# 
# ggplot(vif.df,aes(x=reorder(var,vif),y=vif,fill=group))+
#   geom_bar(stat="identity")+
#   coord_flip()+
#   scale_fill_manual(values=c("#9CBC50","grey60","#5BA79E","#325F92","#E9CF6C","#D55C43"))+
#   theme_classic()+
#   xlab("variable")+
#   geom_hline(yintercept = 10, linetype="dashed", 
#              color = "black", size=0.7)+
#   theme(legend.position="none")
# #scale_fill_manual(values=rev(c("#D66F5D","#D66F5D","grey60","grey60","#5BA79E","#325F92","#E9CF6C","#D55C43","#D55C43","#5BA79E")))
# 
# 
# 
# 
# 
# barplot(vif.cca(nereo.rda))
# 
# plot(nereo.rda, scaling=3)
# 
# #idenitfy rda candidates - snps candidates for local adaptation
# load.rda <- summary(nereo.rda)$species[,1:3]
# #tails (outside normal distribution) are candidates
# hist(load.rda[,1], main="Loadings on RDA1")
# hist(load.rda[,2], main="Loadings on RDA2") 
# hist(load.rda[,3], main="Loadings on RDA3") 
# 
# outliers <- function(x,z){
#   lims <- mean(x) + c(-1, 1) * z * sd(x) ## f.nd loadings +/- z SD from mean loading     
#   x[x < lims[1] | x > lims[2]]           # locus names in these tails
# }
# 
# cand1 <- outliers(load.rda[,1], 3) 
# cand2 <- outliers(load.rda[,2], 3)
# cand3 <- outliers(load.rda[,3], 3)
# 
# nereo.rda.cand <- c(names(cand1), names(cand2), names(cand3)) #names
# length(nereo.rda.cand[duplicated(nereo.rda.cand)]) #no duplicates (on more than 1 axis)
# 
# #closer look at snps
# bgcol  <- ifelse(rdanames %in% nereo.rda.cand, 'gray32', '#00000000')
# snpcol <- ifelse(rdanames %in% nereo.rda.cand, 'red', '#00000000')
# 
# #1 and 2
# plot(nereo.rda, type="n", scaling=3, xlim=c(-0.5,0.5), ylim=c(-0.5,0.5), main="nereo RDA, axes 1 and 2")
# points(nereo.rda, display="species", pch=21, cex=1, col="gray32", bg='#f1eef6', scaling=3)
# points(nereo.rda, display="species", pch=21, cex=1, col=bgcol, bg=snpcol, scaling=3)
# text(nereo.rda, scaling=3, display="bp", col="#0868ac", cex=1)
# 
# #2 and 3
# plot(nereo.rda, type="n", scaling=3, xlim=c(-1,1), ylim=c(-1,1), choices=c(2,3), main="nereo RDA, axes 2 and 3")
# points(nereo.rda, display="species", pch=21, cex=1, col="gray32", bg='#f1eef6', scaling=3, choices=c(2,3))
# points(nereo.rda, display="species", pch=21, cex=1, col=bgcol, bg=snpcol, scaling=3, choices=c(2,3))
# text(nereo.rda, scaling=3, display="bp", col="#0868ac", cex=1, choices=c(2,3))
# 
# #which variables are the most important?
# intersetcor(nereo.rda)[,1:3]
# 
# #               RDA1        RDA2         RDA3
# # SalinityAVG2  0.04229851 -0.842022815 -0.39382484
# # PAR.Summer2   0.56918174 -0.050827540  0.65906136
# # ChlA.Winter   0.14769991  0.046031856 -0.29409376
# # ChlA.Spring  -0.39560984  0.469232487 -0.01815246
# # kd490.Winter  0.46749417  0.295106851 -0.21964747
# # kd490.Spring  0.08787913  0.003939889 -0.21445303
# # lat          -0.95573961 -0.195807751 -0.11037755
# # lon           0.18948244  0.875582403 -0.09351665
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
