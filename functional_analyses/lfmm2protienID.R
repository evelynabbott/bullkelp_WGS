
#Create input for top_go analysis based on SNPs identified by LFMM.R

library(lfmm)
#library(qvalue)
library(stringr)
library(dplyr)
library(vegan)
library(ggplot2)

#load inputs -----
rm(list=ls())
load("pvals4go.RData")
load("qvals4go.RData")

# pvalues<-Nereo.pv.kd490spring.k3$calibrated.pvalue
# qvalues<-Nereo.qv.kd490spring.k3

# pvalues<-Nereo.pv.kd490winter.k3$calibrated.pvalue
# qvalues<-Nereo.qv.kd490winter.k3

pvalues<-Nereo.pv.Parsummer.k3$calibrated.pvalue
qvalues<-Nereo.qv.Parsummer.k3

# pvalues<-Nereo.pv.Salinity.k3$calibrated.pvalue
# qvalues<-Nereo.qv.salinity.k3

# pvalues<-Nereo.pv.WinterSST.k3$calibrated.pvalue
# qvalues<-Nereo.qv.WinterSST.k3


#save these to clear out wd. reload below
save(pvalues,qvalues,file="pqval.Rdata")

#-----------------------------------------------------------------------
rm(list=ls())
load("pqval.Rdata")
sizes.genome<-read.table("sizes.genome",header=F)

# I need to use Scaffold number and snp position which are stored in The 
# pvalue object rownames. I prefer to have that information in a data.frame
# using stringr I will extract that information
#note I have to extract the first row of pvalues, that row has some strange info
SCAFF<-as.numeric(str_extract(rownames(pvalues),regex("(?<=(_))[0-9]+")))
#POS<-as.numeric(str_extract(rownames(pvalues)[-1],regex("(?<=(-))[0-9]+")))
POS=as.numeric(gsub(".*\\.","",rownames(pvalues)))

# pvalues.DF<-data.frame(SCAFF,POS,P.VAL=as.numeric(pvalues)[-1],Q.VAL=qvalues[-1])
pvalues.DF<-data.frame(SCAFF,POS,P.VAL=as.numeric(pvalues),Q.VAL=qvalues)


#I will do the same for the size.genome objects containing the scaffold sizes
#note that the scaffold numbers are larger then the total number of scaffolds
# i.e., the file has 1561 rows and the largest code is 2048

SCAFF.GEN<-as.numeric(str_extract(sizes.genome[,1],regex("(?<=(_))[0-9]+")))
sizes.genome.DF<-data.frame(SCAFF.GEN,SIZE=sizes.genome[,2])

U.scaff.codesPvalue<-unique(pvalues.DF$SCAFF)
#U.scaff.codesPvalue = U.scaff.codesPvalue[1:969] #remove the NA
n.scaffs.in.pval<-length(U.scaff.codesPvalue)

all(U.scaff.codesPvalue%in%sizes.genome.DF$SCAFF.GEN)
#TRUE, just checking that all scaffold exist in the genome

#Make manhattan plots

# for(s in 1:n.scaffs.in.pval){
# 
#   pvalues.DF.temp<-pvalues.DF[pvalues.DF$SCAFF==U.scaff.codesPvalue[s],]
# 
#   LOG10.PVAL<-(-log10(pvalues.DF.temp$P.VAL))
# 
#   XLIM<-sizes.genome.DF$SIZE[sizes.genome.DF$SCAFF.GEN==U.scaff.codesPvalue[s]]
# 
#   plot(x=1:XLIM,ylim=c(0,max(LOG10.PVAL)+2),type="n")
#   SIG<-which(pvalues.DF.temp$Q.VAL < 0.1)
#   notSIG<-which(pvalues.DF.temp$Q.VAL >= 0.1)
#   points(x=pvalues.DF.temp$POS[notSIG],y=LOG10.PVAL[notSIG],col="grey80")
#   points(x=pvalues.DF.temp$POS[SIG],y=LOG10.PVAL[SIG],col="red4")
# 
# 
# }

# for(s in 1:n.scaffs.in.pval){
#   
#   pvalues.DF.temp<-pvalues.DF[pvalues.DF$SCAFF==U.scaff.codesPvalue[s],]
#   
#   LOG10.PVAL<-(-log10(pvalues.DF.temp$P.VAL))
#   
#   XLIM<-sizes.genome.DF$SIZE[sizes.genome.DF$SCAFF.GEN==U.scaff.codesPvalue[s]]
#   
#   plot(x=1:XLIM,ylim=c(0:10),type="n")
#   SIG<-which(pvalues.DF.temp$Q.VAL < 0.1)
#   notSIG<-which(pvalues.DF.temp$Q.VAL >= 0.1)
#   points(x=pvalues.DF.temp$POS[notSIG],y=LOG10.PVAL[notSIG],col="grey80")
#   points(x=pvalues.DF.temp$POS[SIG],y=LOG10.PVAL[SIG],col="red4")
#   
#   
# }


# Write to file the name the scaffold and position of significant snps
#write.table(pvalues.DF[pvalues.DF$Q.VAL<0.1,],"K3_significant_SNPs.txt",sep="\t",quote=F,row.names=F)
write.table(pvalues.DF[pvalues.DF$Q.VAL<0.05,],"K3_significant_SNPs.txt",sep="\t",quote=F,row.names=F)

# Reading the annotation file ----
Nereo.gff<-read.delim("Nerluet1_FilteredModels1_2024-03-15.gff3",header=F,skip=2)

#in the annotation file genes have proteinID field, which is the field that will allow correspondence between annotation file and other databases, e.g., GO terms
#Find rows with the proteinID field
PIDindex<-grep("proteinId",Nereo.gff[,9])

unique(Nereo.gff[PIDindex,3])# we have genes and mRNA with protein IDs

Nereo.gff.PID<-Nereo.gff[PIDindex,]
SCAFF.Nereo.gff.PID<-as.numeric(str_extract(Nereo.gff.PID[,1],regex("(?<=(SCAF_))[0-9]+")))

Nereo.gff.PID<-cbind(SCAFF.Nereo.gff.PID,Nereo.gff.PID)
names(Nereo.gff.PID)[c(1,4,5,6)]<-c("SCAFF","CLASS","START","STOP")

SAMPLES_sig<-read.delim("K3_significant_SNPs.txt")

SCAFF.SNP<-unique(SAMPLES_sig$SCAFF)

#TODO ADD QVALUE TO THE HIT.DF object - match up protein IDs with SNPs
scan.dist<-1000 #can be within 1000 bp of start/stop position

hit.DF<-NULL

for(s in 1:length(SCAFF.SNP)){
  
  tempGFF<-Nereo.gff.PID[Nereo.gff.PID$SCAFF==SCAFF.SNP[s],]
  temp.SNP<-SAMPLES_sig[SAMPLES_sig$SCAFF==SCAFF.SNP[s],]
  
  for(snps in 1:nrow(temp.SNP)){
    d1<-temp.SNP$POS[snps]-tempGFF$START
    d2<-temp.SNP$POS[snps]-tempGFF$STOP
    
    hit.annot.L<-((d1<=0 & abs(d1)<=scan.dist) | (d1>0 & d2<0) | (d2>=0 & d2<=scan.dist))
    
    tempDF<-tempGFF[	hit.annot.L,]
    tempDF<-cbind(rep(temp.SNP$POS[snps],nrow(tempDF)),tempDF)
    names(tempDF)[1]<-"SNP.POS"
    hit.DF<-rbind(hit.DF,tempDF)
    
  }
  print(paste0("Completed Scaffold: ",SCAFF.SNP[s] ))
}

length(unique(hit.DF$SNP.POS))
#Next I run the loop above with scan distance of 1 kbp
length(unique(hit.DF$SNP.POS))

colnames(hit.DF)

hit.DF$proteinId<-as.numeric(str_extract(hit.DF$V9,regex("(?<=(proteinId=))[0-9]+")))

hit.DF.SAMPLES<-hit.DF


# Matching hit protein IDs with information in data base

#Reading the database correspondence files ------------------
#GO TERMS
GO.DF<-read.delim("Nerluet1_FilteredModels1_go_2024-03-15.tab")

#match up protein IDs
GO.DF.SAMPLES<-GO.DF[  GO.DF$proteinId%in%hit.DF.SAMPLES$proteinId, ]
length(unique(GO.DF.SAMPLES$go_term_id))
uniqueGO = unique(GO.DF.SAMPLES$go_name)

#save for top go
#save(hit.DF.SAMPLES,GO.DF.SAMPLES,GO.DF,file="GO.DF.August.Rdata")
#save(hit.DF.SAMPLES,GO.DF.SAMPLES,GO.DF,file="GO.DF.Winter.Rdata")
#save(hit.DF.SAMPLES,GO.DF.SAMPLES,GO.DF,file="GO.DF.Salinity.Rdata")
#save(hit.DF.SAMPLES,GO.DF.SAMPLES,GO.DF,file="GO.DF.Lat.Rdata")
#save(hit.DF.SAMPLES,GO.DF.SAMPLES,GO.DF,file="GO.DF.kd490winter.Rdata")
#save(hit.DF.SAMPLES,GO.DF.SAMPLES,GO.DF,file="GO.DF.PARsummer.Rdata")
#save(hit.DF.SAMPLES,GO.DF.SAMPLES,GO.DF,file="GO.DF.chlAwinter.Rdata")
#save(hit.DF.SAMPLES,GO.DF.SAMPLES,GO.DF,file="GO.DF.kd490spring.Rdata")



