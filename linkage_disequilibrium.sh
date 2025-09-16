



#make sample lists
bcftools query -l WA.Nl.64.F.filtered.vcf.gz > samps
R

samps = read.table("samps")
PopCode<-substr(samps$V1,1,2)

Group1<-c("FW","MC","JT","AH","PT","KR","FB")
Group2<-c("DB","SH","PP","HI","ED","SM","EB")
Group3<-c("LP","PV","SB","DI","SQ")

GenGroupK3<-1:length(samps)
GenGroupK3[PopCode%in%Group1]<-"SJF"
GenGroupK3[PopCode%in%Group2]<-"WB"
GenGroupK3[PopCode%in%Group3]<-"SPS"

sjf.samps = samps[which(substr(samps$V1,1,2) %in% Group1 ), ]
wb.samps = samps[which(substr(samps$V1,1,2) %in% Group2 ), ]
sps.samps = samps[which(substr(samps$V1,1,2) %in% Group3 ), ]

write.table(sjf.samps,file="sjf.samps",col.names = F, row.names = F, quote = F)
write.table(wb.samps,file="wb.samps",col.names = F, row.names = F, quote = F)
write.table(sps.samps,file="sps.samps",col.names = F, row.names = F, quote = F)


bcftools view   --samples-file  sjf.samps    WA.Nl.64.F.filtered.vcf.gz > WA.Nl.64.F.filtered.sjf.vcf.gz #code to subset a VCF file based on a list of file names
bcftools view   --samples-file  wb.samps    WA.Nl.64.F.filtered.vcf.gz > WA.Nl.64.F.filtered.wb.vcf.gz #code to subset a VCF file based on a list of file names
bcftools view   --samples-file  sps.samps    WA.Nl.64.F.filtered.vcf.gz > WA.Nl.64.F.filtered.sps.vcf.gz #code to subset a VCF file based on a list of file names

#index?
#gunzip WA.Nl.64.F.filtered.vcf.gz
# bgzip -c WA.Nl.64.F.filtered.vcf > WA.Nl.64.F.filtered.vcf.gz
# bcftools index WA.Nl.64.F.filtered.vcf.gz

rm *newID*

VCF="WA.Nl.64.F.filtered.sjf.vcf.gz"
#Creating a scaffold map for plink
bcftools view -H $VCF | cut -f 1 | uniq > uniq.txt
#jump to R
R
scaffs<-scan("uniq.txt",what=character())
write.table(data.frame(scaffs,scaffs),file="scaff.map.txt",sep="\t",quote=F,row.names=F,col.names=F)

#Using vcftools to convert vcf to bed, need the scaffold map from above
vcftools --gzvcf $VCF --plink --chrom-map scaff.map.txt --out WA.Nl.64.F.filtered

PED="WA.Nl.64.F.filtered.ped"

#Edit ped to contain population information
less $PED | cut -f 1,2 > updateID.txt

R
ID<-read.table("updateID.txt")
write.table(cbind(ID,substr(ID[,1],1,2),ID[,1]),"IDedit.txt",quote=F,col.names=F,row.names=F)

plink -file WA.Nl.64.F.filtered -aec -update-ids IDedit.txt  -recode -out WA.Nl.64.F.filtered.newID
plink -file WA.Nl.64.F.filtered.newID -aec -update-ids IDedit.txt -make-bed -out WA.Nl.64.F.filtered.newID
#I also repeated the code above for the Inner Puget sound area


BED="WA.Nl.64.F.filtered.newID"
OUTDIR="/scratch/05410/eabbott/kelp/ld.by.pop/sjf/ld_files/send2home"
#Estimating LD
plink -bfile $BED -aec -chr SCAF_1 -r2 -ld-window-kb 1000 -ld-window 100 -ld-window-r2 0 -out $OUTDIR/scaffold1.r0II  
plink -bfile $BED -aec -chr SCAF_2 -r2 -ld-window-kb 1000 -ld-window 100 -ld-window-r2 0 -out $OUTDIR/scaffold2.r0II 
plink -bfile $BED -aec -chr SCAF_3 -r2 -ld-window-kb 1000 -ld-window 100 -ld-window-r2 0 -out $OUTDIR/scaffold3.r0II  
plink -bfile $BED -aec -chr SCAF_4 -r2 -ld-window-kb 1000 -ld-window 100 -ld-window-r2 0 -out $OUTDIR/scaffold4.r0II  
plink -bfile $BED -aec -chr SCAF_5 -r2 -ld-window-kb 1000 -ld-window 100 -ld-window-r2 0 -out $OUTDIR/scaffold5.r0II  
plink -bfile $BED -aec -chr SCAF_6 -r2 -ld-window-kb 1000 -ld-window 100 -ld-window-r2 0 -out $OUTDIR/scaffold6.r0II  
plink -bfile $BED -aec -chr SCAF_7 -r2 -ld-window-kb 1000 -ld-window 100 -ld-window-r2 0 -out $OUTDIR/scaffold7.r0II  
plink -bfile $BED -aec -chr SCAF_8 -r2 -ld-window-kb 1000 -ld-window 100 -ld-window-r2 0 -out $OUTDIR/scaffold8.r0II  
plink -bfile $BED -aec -chr SCAF_9 -r2 -ld-window-kb 1000 -ld-window 100 -ld-window-r2 0 -out $OUTDIR/scaffold9.r0II  
plink -bfile $BED -aec -chr SCAF_10 -r2 -ld-window-kb 1000 -ld-window 100 -ld-window-r2 0 -out $OUTDIR/scaffold10.r0II  
plink -bfile $BED -aec -chr SCAF_11 -r2 -ld-window-kb 1000 -ld-window 100 -ld-window-r2 0 -out $OUTDIR/scaffold11.r0II
plink -bfile $BED -aec -chr SCAF_12 -r2 -ld-window-kb 1000 -ld-window 100 -ld-window-r2 0 -out $OUTDIR/scaffold12.r0II  
plink -bfile $BED -aec -chr SCAF_13 -r2 -ld-window-kb 1000 -ld-window 100 -ld-window-r2 0 -out $OUTDIR/scaffold13.r0II  
plink -bfile $BED -aec -chr SCAF_14 -r2 -ld-window-kb 1000 -ld-window 100 -ld-window-r2 0 -out $OUTDIR/scaffold14.r0II  
plink -bfile $BED -aec -chr SCAF_15 -r2 -ld-window-kb 1000 -ld-window 100 -ld-window-r2 0 -out $OUTDIR/scaffold15.r0II  
plink -bfile $BED -aec -chr SCAF_16 -r2 -ld-window-kb 1000 -ld-window 100 -ld-window-r2 0 -out $OUTDIR/scaffold16.r0II  
plink -bfile $BED -aec -chr SCAF_17 -r2 -ld-window-kb 1000 -ld-window 100 -ld-window-r2 0 -out $OUTDIR/scaffold17.r0II  
plink -bfile $BED -aec -chr SCAF_18 -r2 -ld-window-kb 1000 -ld-window 100 -ld-window-r2 0 -out $OUTDIR/scaffold18.r0II  
plink -bfile $BED -aec -chr SCAF_19 -r2 -ld-window-kb 1000 -ld-window 100 -ld-window-r2 0 -out $OUTDIR/scaffold19.r0II  
plink -bfile $BED -aec -chr SCAF_20 -r2 -ld-window-kb 1000 -ld-window 100 -ld-window-r2 0 -out $OUTDIR/scaffold20.r0II  
plink -bfile $BED -aec -chr SCAF_21 -r2 -ld-window-kb 1000 -ld-window 100 -ld-window-r2 0 -out $OUTDIR/scaffold21.r0II  
plink -bfile $BED -aec -chr SCAF_22 -r2 -ld-window-kb 1000 -ld-window 100 -ld-window-r2 0 -out $OUTDIR/scaffold22.r0II  
plink -bfile $BED -aec -chr SCAF_23 -r2 -ld-window-kb 1000 -ld-window 100 -ld-window-r2 0 -out $OUTDIR/scaffold23.r0II  
plink -bfile $BED -aec -chr SCAF_24 -r2 -ld-window-kb 1000 -ld-window 100 -ld-window-r2 0 -out $OUTDIR/scaffold24.r0II  
plink -bfile $BED -aec -chr SCAF_25 -r2 -ld-window-kb 1000 -ld-window 100 -ld-window-r2 0 -out $OUTDIR/scaffold25.r0II  
plink -bfile $BED -aec -chr SCAF_26 -r2 -ld-window-kb 1000 -ld-window 100 -ld-window-r2 0 -out $OUTDIR/scaffold26.r0II  
plink -bfile $BED -aec -chr SCAF_27 -r2 -ld-window-kb 1000 -ld-window 100 -ld-window-r2 0 -out $OUTDIR/scaffold27.r0II  
plink -bfile $BED -aec -chr SCAF_28 -r2 -ld-window-kb 1000 -ld-window 100 -ld-window-r2 0 -out $OUTDIR/scaffold28.r0II  
plink -bfile $BED -aec -chr SCAF_29 -r2 -ld-window-kb 1000 -ld-window 100 -ld-window-r2 0 -out $OUTDIR/scaffold29.r0II 
plink -bfile $BED -aec -chr SCAF_30 -r2 -ld-window-kb 1000 -ld-window 100 -ld-window-r2 0 -out $OUTDIR/scaffold30.r0II 
plink -bfile $BED -aec -chr SCAF_31 -r2 -ld-window-kb 1000 -ld-window 100 -ld-window-r2 0 -out $OUTDIR/scaffold31.r0II 
plink -bfile $BED -aec -chr SCAF_32 -r2 -ld-window-kb 1000 -ld-window 100 -ld-window-r2 0 -out $OUTDIR/scaffold32.r0II  
plink -bfile $BED -aec -chr SCAF_33 -r2 -ld-window-kb 1000 -ld-window 100 -ld-window-r2 0 -out $OUTDIR/scaffold33.r0II  
plink -bfile $BED -aec -chr SCAF_34 -r2 -ld-window-kb 1000 -ld-window 100 -ld-window-r2 0 -out $OUTDIR/scaffold34.r0II   


#Move to R to plot LD with distance

nano ld_dfs.R 

#write Rscript for TACC -------
# library(tidyverse)
# scaf.stats = list()
# for (i in 1:34){
#   scaf.stats[[i]] = read.table(paste0("scaffold",i,".r0II.ld"),header=T)
#   scaf.stats[[i]]$Dist = scaf.stats[[i]]$BP_B - scaf.stats[[i]]$BP_A
# }

# RII = bind_rows(scaf.stats)
# DistClasses<-cut(RII$Dist,breaks=c(0,100,500,seq(1000,max(RII$Dist),20000)))
# R2classes<-cut(RII$R2,breaks=seq(0,1,0.1))

# calc.means = data.frame("MeanR2" = tapply(RII$R2,DistClasses,mean), "MeanDist" = tapply(RII$Dist,DistClasses,mean))

# save(RII,calc.means,file = "sps_ld.Rdata")

RII_small = RII[sample(nrow(RII), 200000), ]
save(RII_small,calc.means,file="sjf_small.Rdata")
#-------

#run as a job
echo "Rscript ld_dfs.R" > ldrun
ls6_launcher_creator.py -n ldrun -j ldrun -a $allo -e $email -q vm-small -t 00:10:00 -N 1 -w 1
sbatch ldrun.slurm



rsync -v eabbott@ls6.tacc.utexas.edu:/scratch/05410/eabbott/kelp/ld.by.pop/sjf/ld_files/send2home/sjf_small.Rdata .
rsync -v eabbott@ls6.tacc.utexas.edu:/scratch/05410/eabbott/kelp/ld.by.pop/sps/ld_files/send2home/sps_small.Rdata .
rsync -v eabbott@ls6.tacc.utexas.edu:/scratch/05410/eabbott/kelp/ld.by.pop/wb/ld_files/send2home/wb_small.Rdata .







#------------------------------------------------
# Aplot showing LD for all areas together
colors4<-c("orange","tomato","olivedrab",  "lightblue")
s<-1



#R script
library(tidyverse)
for(s in 1:34){
RII.outer<-read.table(paste0("scaffold",s,".seattle.r0II.ld"),header=T)
RII.outer$Dist<-RII.outer$BP_B-RII.outer$BP_A
DistClasses.outer<-cut(RII.outer$Dist,breaks=c(0,100,500,seq(1000,max(RII.outer$Dist),10000)))
R2classes.outer<-cut(RII.outer$R2,breaks=seq(0,1,0.1))
MeanR2.outer<-tapply(RII.outer$R2,DistClasses.outer,mean,na.rm=T)
MeanDist.outer<-tapply(RII.outer$Dist,DistClasses.outer,mean,na.rm=T)
}

save(list = ls(pattern = 'outer'), file = "outer_res.RData")
##

echo "Rscript outer.R" > ldrun
ls6_launcher_creator.py -n ldrun -j ldrun -a $allo -e $email -q vm-small -t 00:10:00 -N 1 -w 1
sbatch ldrun.slurm



###

png(file="file.png",width=1750, height=1500)
plot(MeanDist.outer,MeanR2.outer,xlim=c(0,120000))
points(RII.outer$Dist,RII.outer$R2,col="grey85",cex=0.7)


dev.off()

#rsync sjf
scp eabbott@ls6.tacc.utexas.edu:/scratch/05410/eabbott/kelp/ld.by.pop/sjf/ld_files/send2home/ .

scp -v eabbott@ls6.tacc.utexas.edu:/scratch/05410/eabbott/kelp/ld.by.pop/sps/ld_files/ .











