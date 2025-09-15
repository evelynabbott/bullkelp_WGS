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

load("lfmm64_input.Rdata")

pairs(Env.dataP)
Env.dataP = Env.dataP[,-c(4,7,8,9)]

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

Nereo.lfmmPC1.k3 <- lfmm_ridge(Y=gen.ImpN.df, X=pred.PC1, K=3) ## c.ange K as you see fit
Nereo.lfmmPC2.k3 <- lfmm_ridge(Y=gen.ImpN.df, X=pred.PC2, K=3) ## c.ange K as you see fit

#LFMM per variable (only those identified by RDAforest)
Nereo.lfmmSST.Salinity.k3 <- lfmm_ridge(Y=gen.ImpN.df, X=Env.dataP$SalinityAVG2, K=3)
Nereo.lfmmSST.Winter.k3 <- lfmm_ridge(Y=gen.ImpN.df, X=Env.dataP$SST.month.Winter, K=3)
Nereo.lfmm.Parsummer.k3 <- lfmm_ridge(Y=gen.ImpN.df, X=Env.dataP$PAR.Summer2, K=3)
Nereo.lfmm.kd490winter.k3 <- lfmm_ridge(Y=gen.ImpN.df, X=Env.dataP$kd490.Winter, K=3)
Nereo.lfmm.kd490spring.k3 <- lfmm_ridge(Y=gen.ImpN.df, X=Env.dataP$kd490.Spring, K=3)

#PCs pvalues
Nereo.pv.pc1.k3 <- lfmm_test(Y=gen.ImpN.df, X=pred.PC1, lfmm=Nereo.lfmmPC1.k3, calibrate="gif")
Nereo.pv.pc2.k3 <- lfmm_test(Y=gen.ImpN.df, X=pred.PC2, lfmm=Nereo.lfmmPC2.k3, calibrate="gif")

#variable pvals
Nereo.pv.Salinity.k3<- lfmm_test(Y=gen.ImpN.df, X=Env.dataP$SalinityAVG, lfmm=Nereo.lfmmSST.Salinity.k3, calibrate="gif")
Nereo.pv.WinterSST.k3<- lfmm_test(Y=gen.ImpN.df, X=Env.dataP$SST.month.Winter, lfmm=Nereo.lfmmSST.Winter.k3, calibrate="median+MAD")
Nereo.pv.Parsummer.k3<- lfmm_test(Y=gen.ImpN.df, X=Env.dataP$PAR.Summer2, lfmm=Nereo.lfmm.Parsummer.k3, calibrate="gif")
Nereo.pv.kd490winter.k3<- lfmm_test(Y=gen.ImpN.df, X=Env.dataP$kd490.Winter, lfmm=Nereo.lfmm.kd490winter.k3, calibrate="gif")
Nereo.pv.kd490spring.k3<- lfmm_test(Y=gen.ImpN.df, X=Env.dataP$kd490.Spring, lfmm=Nereo.lfmm.kd490spring.k3, calibrate="gif")

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


# 
#q vals ------
#K3
Nereo.qv.pc1.k3 <- qvalue(Nereo.pv.pc1.k3$calibrated.pvalue)$qvalues
Nereo.qv.pc2.k3 <- qvalue(Nereo.pv.pc2.k3$calibrated.pvalue)$qvalues

Nereo.qv.salinity.k3 <- qvalue(Nereo.pv.Salinity.k3$calibrated.pvalue)$qvalues
Nereo.qv.WinterSST.k3 <- qvalue(Nereo.pv.WinterSST.k3$calibrated.pvalue)$qvalues
Nereo.qv.Parsummer.k3 <- qvalue(Nereo.pv.Parsummer.k3$calibrated.pvalue)$qvalues
Nereo.qv.kd490winter.k3 <- qvalue(Nereo.pv.kd490winter.k3$calibrated.pvalue)$qvalues
Nereo.qv.kd490spring.k3 <- qvalue(Nereo.pv.kd490spring.k3$calibrated.pvalue)$qvalues

#save these files to clear up space in env
save(list = ls(pattern = '.pv.'), file = "pvals4go.RData")
save(list = ls(pattern = '.qv.'), file = "qvals4go.RData")


#These can be used to extract the loci under section associated with different 
#PC1
## how many SNPs have an FDR < 10%?
length(which(Nereo.qv.pc1.k3 < 0.1)) # 1065
length(which(Nereo.qv.pc2.k3 < 0.1)) # 221

length(which(Nereo.qv.WinterSST.k3 < 0.05)) # 12249
length(which(Nereo.qv.salinity.k3 < 0.05)) # 141 *** 1
length(which(Nereo.qv.kd490winter.k3 < 0.05)) # 413
length(which(Nereo.qv.kd490spring.k3 < 0.05)) # 48 * 3
length(which(Nereo.qv.Parsummer.k3 < 0.05)) # 1725 ** 2

