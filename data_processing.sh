
#creating a file with all file names
# in R
NamesTrim<-grep("R1_001",file.names,value=T)
NamesTrim<-sub("R1_001.fastq.gz","",NamesTrim)

write(NamesTrim,"NamesTrim.txt")


# Jump out of R
# conda activate bioinfo ( I use conda environments)

cat NamesTrim.txt | parallel -j 18 fastqc {}R1_001.fastq.gz {}R2_001.fastq.gz


#used instruction from biostar book
#installed fastp and used the code below to find the adapter
#Note that only works in a single end for pair end sequencing (we are only detecting the adapter"
#We found the adapter below (there is more to the output) now we can create a fasta file with the adapter

fastp -i AH10MD6_S28_L004_R1_001.fastq.gz AH10MD6_S28_L004_R1_001.trimmed.fq

echo ">illumina" > adapter.fa
echo "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA" >>adapter.fa

fastp --cut_tail --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
       -i AH10MD6_S28_L004_R1_001.fastq.gz -o AH10MD6_S28_L004_R1_001.trim.fq.gz \
       -I AH10MD6_S28_L004_R2_001.fastq.gz -O AH10MD6_S28_L004_R2_001.trim.fq.gz


# Doing all in parallel
#Create the appropriate file with sample names
#back in R
file.names1<-sub("R1_001.fastq.gz","",file.names)
file.names2<-sub("R2_001.fastq.gz","",file.names1)
write(grep("L00",file.names2,value=T),"NamesTrim2.txt")  #NamesTrim2.txt contains the sample names ready for parallel

cat NamesTrim2.txt | parallel -j 18 fastp --cut_tail --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
       -i {}R1_001.fastq.gz -o {}R1_001.trim.fq.gz \
       -I {}R2_001.fastq.gz -O {}R2_001.trim.fq.gz


# Reference genome
# indexing the reference genome
#For JGI (annotated) version
bwa index Nerluet1_AssemblyScaffolds_2024-03-15.fasta

#local variable
REF=/media/Backup6t/WA_WGS_Nereo/reference/JGI/Nerluet1_AssemblyScaffolds_2024-03-15.fasta


#NamesTrim.txt contains the sample names without the common sufixes
cat NamesTrim.txt | parallel -j 18 bwa mem -M -t 1 $REF {}R1_001.trim.fq.gz {}R2_001.trim.fq.gz \> ./align/{}.sam 


#Let's get flagstats stats for sam sequences
R
samFiles<-list.files()
#first file is the folder nae flagstat
write(samFiles[-1],"samFiles.txt")

cat samFiles.txt | parallel -j 12 samtools flagstat {} \> ./flagstat/{}.flagstat

#Compressing to bam format
cat SamFiles | parallel -j 18 samtools view {}_.sam -b -o {}.bam

#sorting bam
cat SamFiles | parallel -j 18 samtools sort ./bam/{}.bam -o {}_sort.bam


#Next step if to merge the sorted bam biles from the differente sequence lanes (3 sequences lanes wre used 4, 5, and 6)

#but first I need to create a list of files and check that the files that I have on disk all have three file per samples
# I did that in R and commented out the code here for simplicity

#move sorted.bam files into new folder "sorted.bam"
#R
#sorted.bam.files<-list.files()
#sorted.bam.files<-sorted.bam.files[grep(sorted.bam.files,pattern="_sort.bam")]

#Mtemp<-matrix(
#unlist(strsplit(sorted.bam.files,split="_")),
#byrow=T,ncol=4)

#FileNames<-paste(Mtemp[,1],Mtemp[,2],sep="_")

#all(table(FileNames)==3) #This showed that all 34 files have the 3 sequencing lane files
#get the unique file names
#write(unique(FileNames),"mergeFiles.txt")
#back to the shell

#merge the three files into one
cat mergeFiles.txt | parallel -j 20 samtools merge {}.bam {}_L004_sort.bam {}_L005_sort.bam {}_L006_sort.bam

#Don't forget to sort again after merging
cat mergeFiles.txt | parallel -j 20 samtools sort {}.bam -o {}_sort.bam



# Finally we need to remove duplicate reads with picard tools #the j flag sets how many files parallel will do at a time

cat mergeFiles.txt | parallel --verbose -j 20 \
java -Xmx1g -jar /home/filipe/miniconda3/envs/ddrad/share/picard-2.25.6-0/picard.jar \
 MarkDuplicates REMOVE_DUPLICATES=true \
 ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT \
 MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 \
 INPUT={}_sort.bam \
 OUTPUT={}_sort.rmd.bam \
 METRICS_FILE={}_sort.rmd.bam.metrics


# Now we need to index all bam files again 
cat mergeFiles.txt | parallel -j 21 samtools index {}_sort.rmd.bam


#Calling variatnts
indexing the genome again


samtools faidx ./reference/JGI/Nerluet1_AssemblyScaffolds_2024-03-15.fasta

REF=./reference/Nerluet1_AssemblyScaffolds_2024-03-15.fasta

#Calling variants, this is the slow part. Note that I am calling haploids, which should be faster.
bcftools mpileup -a AD,DP,SP -O u -f $REF *_sort.rmd.bam | bcftools call --ploidy 1 -f GQ,GP -m -O z -v -o ./WA.Nl.100.JGIannoted.vcf.gz


#Analzyizing the vcf 
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' WA.Nl.45.vcf.gz | head
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%TGT]\n'  WA.Nl.100.JGIannoted.vcf.gz | head

#How many calls do we have
bcftools view -H   WA.Nl.100.JGIannoted.vcf.gz | wc -l

VCF=WA.Nl.100.JGIannoted.vcf.gz

OUT_F=./vcftools/WA.Nl.100.vcf

vcftools --gzvcf $VCF --freq2 --out $OUT_F --max-alleles 2
vcftools --gzvcf $VCF --depth --out $OUT_F
vcftools --gzvcf $VCF --site-mean-depth --out $OUT_F
vcftools --gzvcf $VCF --site-quality --out $OUT_F
vcftools --gzvcf $VCF --missing-indv --out $OUT_F
vcftools --gzvcf $VCF --het --out $OUT_F  #this one is realy not relevant for haploids

# Moved to R to start analyzing the output of vcftools
# Variants missed per sample
imiss.Full<-read.delim("./vcftools/WA.Nl.100.vcf.imiss")


summary(imiss.Full$F_MISS)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.03336 0.04676 0.04902 0.05057 0.05278 0.07472   #for  45 samples

#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.03581 0.04971 0.05236 0.06105 0.05776 0.44275  #for full 100 samples

#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.03561 0.04958 0.05278 0.06296 0.05926 0.44359 #for full 100 samples new mapping to the annotated genome


lqual.Full<-read.delim("./vcftools/WA.Nl.100.vcf.lqual")

summary(lqual.Full)

 # CHROM                POS                QUAL       
 # Length:5758391     Min.   :       5   Min.   :  3.01  #for  45 samples
 # Class :character   1st Qu.: 2276638   1st Qu.:209.00  
 # Mode  :character   Median : 5151044   Median :578.00  
 #                    Mean   : 5825108   Mean   :579.82  
 #                    3rd Qu.: 8569722   3rd Qu.:999.00  
 #                    Max.   :21866497   Max.   :999.00 
                    
 #          CHROM                POS                QUAL       
 # Length:7176418     Min.   :       5   Min.   :  3.01  #for full 100 samples
 # Class :character   1st Qu.: 2311377   1st Qu.:208.00  
 # Mode  :character   Median : 5175876   Median :603.00  
 #                    Mean   : 5846724   Mean   :585.36  
 #                    3rd Qu.: 8585934   3rd Qu.:999.00  
 #                    Max.   :21866510   Max.   :999.00           
                    
  hist(  lqual.Full$QUAL    )  
  # Doing really good here
  rm(lqual.Full)  # removing the 620Mb object
            
 ldepth.F<-read.delim("./vcftools/WA.Nl.100.vcf.ldepth.mean")
hist( ldepth.F$MEAN_DEPTH)
# summary( ldepth.F$MEAN_DEPTH)
#    summary( ldepth.F$MEAN_DEPTH)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#   0.02222   9.13333  11.91110  11.12018  13.31110 252.93300   #for full 100 samples #low coverage

# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    0.01    8.60   11.10   10.37   12.34  276.80 #for full 100 samples new mapping to the annotated genome



hist( ldepth.F$MEAN_DEPTH,xlim=c(0,100),breaks=seq( 0,max( ldepth.F$MEAN_DEPTH),1))

hist( ldepth.F$MEAN_DEPTH,xlim=c(0,100),breaks=277)
dev.copy2pdf(file="./vcftools/plots/ldepth mean.pdf")

ind_depth <- read.delim("./vcftools/WA.Nl.100.vcf.idepth", 
                        col.names = c("ind", "nsites", "depth"), skip = 1)

#This file is good to filter individuals with poor coverage
ind_depth$ind[ind_depth$depth<10]  #35 samples have fewer than 10X cover
 

#counts per pop for inds with less than 10X

table(substr(ind_depth$ind[ind_depth$depth<10],1,2))
# H DI ED FB FW JT KR MC SB SQ 
#  3  7  1  1  3  4  1  3  9  2 
 
# AH DI ED FB FW JT KR MC SB SQ 
#  3  7  1  1  3  4  1  3 10  2  
 
 #counts per pop for inds with more or equal than 10X

table(substr(ind_depth$ind[ind_depth$depth>=10],1,2)) #samples with 10X or more coverage
# AH DB DI EB ED FB FW HI JT KR LP MC PP PT PV SB SH SM SQ 
#  3  5  4  5  4  3  1  2  2  3  5  2  3  3  4  3  5  5  2 
#write to file individuals to remove from vcf file

write(ind_depth$ind[ind_depth$depth<10],"./vcftools/samples2remove")

#Minor allele frequency
frq.Full<-read.delim("./vcftools/WA.Nl.100.vcf.frq",
		 col.names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"),
		  skip = 1)

maf<-apply(frq.Full[,c("a1","a2")],1,min)

length(maf)
[1] 6 949 457 total number of calls

#For 45 samples
 sum(maf<0.05)
# [1]2 014 323    #45 samples, 2M out of 5.6M that have minor allele freq below 5%
# [1] 3 140 223    #100 samples, 3.1M SNPs at 5% MAF
# [1] 3 138 963    #for full 100 samples new mapping to the annotated genome



sum(maf<0.1)
# [1] 2 709 962  #increases to 2.7M if we increase maf filtering to >10%
# [1] 4 050 347   #100 samples, 3.1M SNPs at 5% MAF
#  4 047 932    #for full 100 samples new mapping to the annotated genome
 


#Back on the shell  Lety's remove the individuals with low coverage
#First remove individuals with low sequencing depth (the list created above above)
bcftools view -S ^./vcftools/samples2remove  WA.Nl.100.JGIannoted.vcf.gz > WA.Nl.64.F.vcf
#compress again
bcftools view -O z -o WA.Nl.64.F.vcf.gz WA.Nl.64.F.vcf
# index vcf
bcftools index WA.Nl.64.F.vcf.gz

VCF_IN=WA.Nl.64.F.vcf.gz
VCF_OUT=WA.Nl.64.F.filtered.vcf.gz

#applying filters

MAF=0.045
MISS=0.9
QUAL=30  
MIN_DEPTH=10
MAX_DEPTH=40


vcftools --gzvcf $VCF_IN --remove-indels --maf $MAF --max-missing $MISS --minQ $QUAL --min-meanDP $MIN_DEPTH --max-meanDP $MAX_DEPTH --minDP $MIN_DEPTH --maxDP $MAX_DEPTH --recode --stdout | gzip -c > $VCF_OUT

#listing sample names is a vcf file
bcftools query -l WA.Nl.64.F.filtered.vcf.gz


####################################################
############ ADMIXTURE WITH NGSADMIX ###############
####################################################

ls *bam > bams

#Identity by state (IBS) -------------------------------------
export GRate=0.1 
export MI=46 #aim for ~75% of samples. lower MI if the output includes less than 10,000 sites
FILTERS0='-minInd $MI -uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 20 -snp_pval 1e-5 -dosnpstat 1 -doHWE 1 -hetbias_pval 1e-3 -skipTriallelic 1'
TODO0='-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doPost 1 -doGlf 2'
echo 'export NIND=`cat bams | wc -l`; export MI=`echo "($NIND*$GRate+0.5)/1" | bc`' >calc1
echo "source calc1 && angsd -b bams -GL 1 $FILTERS0 $TODO0 -P 12 -out myresult && Rscript ~/bin/detect_clones.R bams myresult.ibsMat 0.04">a1
ls6_launcher_creator.py -j a1 -n a1 -a $allo -e $email -t 10:00:00 -w 1 -q gpu-a100
sbatch a1.slurm

scp -v eabbott@ls6.tacc.utexas.edu:/scratch/05410/eabbott/kelp/bams_mort/hctree.pdf .

#Run ngs_admix ---------------------------------
#`seq 2 6` is the range of potential population numbers
>ngs
for K in `seq 2 6`
do echo "NGSadmix -likes myresult.beagle.gz -K $K -P 10 -o pops_k${K}" >> ngs
done

ls6_launcher_creator.py -j ngs -n ngs -a $allo -e $email -t 2:00:00 -N 2 -w 2 -q development
sbatch ngs.slurm

#scp to local for figure plotting in R
scp -v eabbott@ls6.tacc.utexas.edu:/scratch/05410/eabbott/kelp/bams_mort/bams.nr .
scp -v eabbott@ls6.tacc.utexas.edu:/scratch/05410/eabbott/kelp/bams_mort/myresult.ibsMat .
scp -v eabbott@ls6.tacc.utexas.edu:'/scratch/05410/eabbott/kelp/bams_mort/*qopt' .




