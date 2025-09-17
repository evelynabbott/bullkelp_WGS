#SnpEff - match variants with protein IDs

library(tidyverse)
library(ggplot2)
library(gridExtra)
library(ggpubr)

#snp function
#let's match up the SNPEff results with the annotation files
rm(list=ls())
load("impact_snpbyind.Rdata") #snp impact
sizes.genome<-read.table("sizes.genome",header=F)

impact.all = impact.all[impact.all$impact != "MODIFIER",]

SCAFF<-as.numeric(gsub("\\..*","",impact.all$pos))
POS=as.numeric(gsub(".*\\.","",impact.all$pos))

pvalues.DF<-data.frame(SCAFF,POS)

SCAFF.GEN<-as.numeric(str_extract(sizes.genome[,1],regex("(?<=(_))[0-9]+")))
sizes.genome.DF<-data.frame(SCAFF.GEN,SIZE=sizes.genome[,2])

U.scaff.codesPvalue<-unique(pvalues.DF$SCAFF)
n.scaffs.in.pval<-length(U.scaff.codesPvalue)

all(U.scaff.codesPvalue%in%sizes.genome.DF$SCAFF.GEN)

Nereo.gff<-read.delim("Nerluet1_FilteredModels1_2024-03-15.gff3",header=F,skip=2) #annotation gff

#in the annotation file genes have proteinID field, which is the field that will allow correspondence between annotation file and other databases, e.g., GO terms
#Find rows with the proteinID field
PIDindex<-grep("proteinId",Nereo.gff[,9])

unique(Nereo.gff[PIDindex,3])# we have genes and mRNA with protein IDs

Nereo.gff.PID<-Nereo.gff[PIDindex,]
SCAFF.Nereo.gff.PID<-as.numeric(str_extract(Nereo.gff.PID[,1],regex("(?<=(SCAF_))[0-9]+")))

Nereo.gff.PID<-cbind(SCAFF.Nereo.gff.PID,Nereo.gff.PID)
names(Nereo.gff.PID)[c(1,4,5,6)]<-c("SCAFF","CLASS","START","STOP")

SAMPLES_sig<-pvalues.DF

SCAFF.SNP<-unique(SAMPLES_sig$SCAFF)

#compare variant locations with protein ID start/stop sites
scan.dist<-0 #must be exact match, i.e. variant must fall within coding region, no flanking regions

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
#851
length(unique(hit.DF$SNP.POS))

colnames(hit.DF)

hit.DF$proteinId<-as.numeric(str_extract(hit.DF$V9,regex("(?<=(proteinId=))[0-9]+")))
impact.prot = data.frame("pos" = paste0(hit.DF$SCAFF,".",hit.DF$SNP.POS),"proteinId" = hit.DF$proteinId) %>% distinct()
impact.all = full_join(impact.prot,impact.all,by = "pos") #dataframe including impact, location, protein ID, pop, cluster and individual
impact.prot = impact.all

save(impact.prot,file="impact.prot.Rdata")

#how many high impact coding variants?
# nrow(impact.high[which(impact.high$region=="PS"),])
# nrow(impact.high[which(impact.high$region=="WI"),])
# nrow(impact.high[which(impact.high$region=="SJF"),])


