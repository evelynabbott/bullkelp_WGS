

#make dataframes for plotting lfmm rda & identifying SNPs

# if(!requireNamespace("qvalue", quietly = TRUE)) {  
# 	if (!requireNamespace("BiocManager", quietly = TRUE))
# 		install.packages("BiocManager")
# 	BiocManager::install(version = "3.14")
# 	BiocManager::install("qvalue")
# }
# if(!requireNamespace("lfmm", quietly = TRUE)) {  
# 	remotes::install_github("bcm-uga/lfmm")
# }
# 
# devtools::install_github("gavinsimpson/ggvegan")


library(devtools)
library(ggvegan)
library(lfmm)
library(qvalue)
library(LandGenCourse)
library(stringr)
library(plyr)
library(vegan)
library(ggplot2)


echo "$HOME/bin/bcftools-1.3.1/bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' WA.Nl.64.F.filtered.vcf.gz > WA.Nl.64.F.filtered.vcf.1.0" > make_table
ls6_launcher_creator.py -j make_table -n make_table -a $allo -e $email -t 00:30:00 -N 1 -w 1 -q vm-small
sbatch make_table.slurm

#make sample names table
$HOME/bin/bcftools-1.3.1/bcftools query -l WA.Nl.64.F.filtered.vcf.gz > sample.names

R
rm(list=ls())
genos.1.0<-read.table("WA.Nl.64.F.filtered.vcf.1.0") # All Markers
table(nchar(genos.1.0[,4]))
genos.1.0<-genos.1.0[nchar(genos.1.0[,4])==1,]
genos.1.0G<-genos.1.0[,5:ncol(genos.1.0)]
genos.1.0G[genos.1.0G=="./."]<-NA
genos.1.0G<-t(genos.1.0G)
colnames(genos.1.0G)<-paste(genos.1.0[,1],genos.1.0[,2],sep="-")
length(rownames(genos.1.0G)) #83
save(genos.1.0G,file="genos.1.0G.64.Rdata")

colnames(genos.1.0G) #looks good
head(rownames(genos.1.0G),5)

#scp eabbott@ls6.tacc.utexas.edu:/scratch/05410/eabbott/kelp/lfmm64/genos.1.0G.64.Rdata .

#----------------------------------------------
rm(list=ls())
load("genos.1.0G.64.Rdata")
snpnames = colnames(genos.1.0G)
#save(snpnames,file="snpnames.Rdata")

#extract site names from bam list
#load("env.64")
Sample.names=read.table("sample.names")
Sample.Names = as.data.frame(Sample.names)
PopCode<-substr(Sample.names$V1,1,2)
table(PopCode)
# AH DB DI EB ED FB FW HI JT KR LP MC PP PT PV SB SH SM SQ 
# 4  5  4  5  4  3  1  2  2  3  5  2  3  3  4  3  5  5  2

#Based on admixture see plot in admixture/plots/admicture WA Nereo 65 samples K3.pdf
Group1<-c("FW","MC","JT","AH","KR")
Group2<-c("DB","SH","PP","HI","ED","SM","EB","PT","FB")
Group3<-c("LP","PV","SB","DI","SQ")

#The Vector GenGroupK3 will have the membership code for each site to one of the three genetic groups
GenGroupK3<-1:length(Sample.Names)
GenGroupK3[PopCode%in%Group1]<-1
GenGroupK3[PopCode%in%Group2]<-2
GenGroupK3[PopCode%in%Group3]<-3

#Now we can input again but based on gen group membership
sum(is.na(genos.1.0G))
#so that I can order again as originally ordered
genos.1.0G<-cbind(1:nrow(genos.1.0G),genos.1.0G)

genos.1.0G.K1<-genos.1.0G[GenGroupK3==1,]
genos.1.0G.K2<-genos.1.0G[GenGroupK3==2,]
genos.1.0G.K3<-genos.1.0G[GenGroupK3==3,]


gen.impK1 <- apply(genos.1.0G.K1, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
#sum(is.na(genos.1.0G.K1))
#sum(is.na(gen.impK1)) # 
save(gen.impK1,file="gen.impK1")

gen.impK2 <- apply(genos.1.0G.K2, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
sum(is.na(genos.1.0G.K2)) # 
sum(is.na(gen.impK2)) # 
save(gen.impK2,file="gen.impK2")

gen.impK3 <- apply(genos.1.0G.K3, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
sum(is.na(genos.1.0G.K3)) #
sum(is.na(gen.impK3)) #
save(gen.impK3,file="gen.impK3")

dim(gen.impK3)

#rebind and reorder
gen.imp<-rbind(gen.impK1,gen.impK2,gen.impK3)
gen.imp<-gen.imp[order(as.numeric(gen.imp[,1])),] #reordered

#fix names
fix.names = gsub("\\_.*","",Sample.names$V1)
rownames(gen.imp) = fix.names
gen.imp=gen.imp[,-c(1)]

save(gen.imp,file="gen.imp64.Rdata")

dim(gen.imp) #

gen.impN<-apply(gen.imp,2,as.numeric)
gen.ImpN.df<-data.frame(gen.impN)
#names(gen.ImpN.df)<-colnames(genos.1.0G)
gen.pca <- rda(gen.ImpN.df, scale=T)

dim(gen.ImpN.df)

save(gen.ImpN.df,gen.pca,Sample.Names,file="genImp64.Rdata")




