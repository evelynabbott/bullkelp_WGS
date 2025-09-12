#################### Using TopGO ----

#BiocManager::install("topGO")
library(topGO)
library(gridExtra)

# Following a different tutorial (https://gist.github.com/slavailn/dcff753cc32fd9e053590308c831b057)
rm(list=ls())
input="winter"
#load("GO.DF.kd490winter.Rdata")
#load("GO.DF.PARsummer.Rdata")
#load("GO.DF.Salinity.Rdata")
#load("GO.DF.kd490spring.Rdata")
load("GO.DF.Winter.Rdata")

#save(hit.DF.SAMPLES,file="hit.DF.Rdata")

#First set create a gene2Go table


PID.u<-unique(GO.DF$proteinId)
Genes2GO.M<-matrix(ncol=2,nrow=length(PID.u))

for(p in 1:length(PID.u)){
  Genes2GO.M[p,1]<-as.character(PID.u[p])
  Genes2GO.M[p,2]<-paste(GO.DF$go_acc[GO.DF$proteinId==PID.u[p]],collapse=",")
}

#Note that I am using as gene names protein-ID
write.table(Genes2GO.M,"Genes2GO.M.txt",row.names=F,quote=F,sep="\t") #the remove header!
geneID2GO<-readMappings("Genes2GO.M.txt")
geneNames<-names(geneID2GO)


# Get the list of genes of interest
myInterestingGenes <- unique(GO.DF.SAMPLES$proteinId)   #Salinity

geneList <- factor(as.integer(geneNames %in% myInterestingGenes))   #PC1.K4
names(geneList) <- geneNames

#head(geneList)

#Create geneList above for the three different sets (The geneList object is recreated above for each set)
GOdata.BP <- new("topGOdata", ontology = "BP", allGenes = geneList,
                 annot = annFUN.gene2GO, gene2GO = geneID2GO)

GOdata.MF <- new("topGOdata", ontology = "MF", allGenes = geneList,
                 annot = annFUN.gene2GO, gene2GO = geneID2GO)

GOdata.CC <- new("topGOdata", ontology = "CC", allGenes = geneList,
                 annot = annFUN.gene2GO, gene2GO = geneID2GO)



resultTopGO.elim.MF <- runTest(GOdata.MF, algorithm = "elim", statistic = "Fisher" )
resultTopGO.elim.CC <- runTest(GOdata.CC, algorithm = "elim", statistic = "Fisher" )
resultTopGO.elim.BP <- runTest(GOdata.BP, algorithm = "elim", statistic = "Fisher" )


#PC1.K4

#any
allRes.BP <- GenTable(GOdata.BP, elimKS = resultTopGO.elim.BP,
                              orderBy = "elimKS", 
                              topNodes = 356)

allRes.MF <- GenTable(GOdata.MF, elimKS = resultTopGO.elim.MF,
                              orderBy = "elimKS", 
                              topNodes = 356)

allRes.CC <- GenTable(GOdata.CC, elimKS = resultTopGO.elim.CC,
                              orderBy = "elimKS",
                              topNodes = 303)



#Writing topGO result tables to file
#PC1.K4
#write.table(allRes, file = paste0("topGO_results_",input,".all.txt"), sep = "\t", quote = F, col.names = T, row.names = F)

write.table(allRes.BP, file = paste0("topGO_results_",input,".BP.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
write.table(allRes.MF, file = paste0("topGO_results_",input,".MF.txt"), sep = "\t", quote = F, col.names = T, row.names = F)
write.table(allRes.CC, file = paste0("topGO_results_",input,".CC.txt"), sep = "\t", quote = F, col.names = T, row.names = F)


# Bar plots ----------

#pval<-as.numeric(allRes$elimKS)
pval<-as.numeric(allRes.BP$elimKS)
par(mar=c(5.1 ,22.1, 4.1 ,2.1))
barplot(rev(-log10(pval)[pval<0.01]),horiz=T,las=1,main="Biological Process",xlab="Enrichment",
				names.arg=rev(allRes.BP$Term[pval<0.01]))

pval<-as.numeric(allRes.MF$elimKS)
par(mar=c(5.1 ,22.1, 4.1 ,2.1))
barplot(rev(-log10(pval)[pval<0.01]),horiz=T,las=1,main="Molecular Function",xlab="Enrichment",
				names.arg=rev(allRes.MF$Term[pval<0.01]))

pval<-as.numeric(allRes.CC$elimKS)
par(mar=c(5.1 ,22.1, 4.1 ,2.1))
barplot(rev(-log10(pval)[pval<0.01]),horiz=T,las=1,main="Cellular Component",xlab="Enrichment",
				names.arg=rev(allRes.CC$Term[pval<0.01]))



#subset GO terms which are enriched p < 0.01
bp.sig=allRes.BP[which(allRes.BP$elimKS < 0.05),]
bp.sig$cat = "bp"
cc.sig=allRes.CC[which(allRes.CC$elimKS < 0.05),]
cc.sig$cat = "cc"
mf.sig=allRes.MF[which(allRes.MF$elimKS < 0.05),]
mf.sig$cat = "mf"
sig.gos=rbind(bp.sig,cc.sig,mf.sig)

length(allRes.BP[which(allRes.BP$elimKS < 0.05),])
length(allRes.CC[which(allRes.CC$elimKS < 0.05),])
length(allRes.MF[which(allRes.MF$elimKS < 0.05),])


#save for downstream
assign(paste0(input,"_sig.gos"),sig.gos)
#save(winter_sig.gos,file="winter_sig.gos.cats.Rdata")
#save(kd490winter_sig.gos,file="kd490winter_sig.gos.cats.Rdata")
#save(parsummer_sig.gos,file="PARsummer_sig.gos.cats.Rdata")
#save(salinity_sig.gos,file="salinity_sig.gos.cats.Rdata")
#save(kd490spring_sig.gos,file="kd490spring_sig.gos.cats.Rdata")


#barplots ---------------
rm(list=ls())
load("salinity_sig.gos.cats.Rdata")

df.plot.s = salinity_sig.gos[which(salinity_sig.gos$elimKS <= 0.01),]
# df.plot.s$enrichment = rev(-log10(as.numeric(df.plot.s$elimKS)))
df.plot.s$enrichment = -log10(as.numeric(df.plot.s$elimKS))
colnames(df.plot.s)[2] <- "term"
colnames(df.plot.s)[7] <- "process"

#add separator bars
df.plot.s[6,] = c("","","","","","","none1",1)
#df.plot.s[7,] = c("","","","","","","none2",1)
df.plot.s$enrichment = as.numeric(df.plot.s$enrichment)


df.plot.s$procnum = ifelse(df.plot.s$process == "bp",1,
                      ifelse(df.plot.s$process == "none1",2,
                             ifelse(df.plot.s$process == "cc",3,
                                    ifelse(df.plot.s$process == "none2",4,5))))

a=df.plot.s %>%
  mutate(reord = as.numeric(procnum) + as.numeric(enrichment),
         term = fct_reorder(term, reord, .desc = F)) %>% 
  ggplot(aes(y =  reorder(term, -procnum),  x = enrichment, fill = process)) +
  geom_bar(stat = "identity",  width = 0.8,alpha=0.75)+
  scale_fill_manual(values = c("bp" = "#D59935", "mf" = "#BE3B38"),na.value = NA)+
  theme_minimal()+
  ylab("")+
  xlab("")+
  ggtitle("salinity")+
  theme(legend.position = "none")+
  xlim(0,3.3)

#PAR
load("PARsummer_sig.gos.cats.Rdata")

df.plot.p = parsummer_sig.gos[which(parsummer_sig.gos$elimKS <= 0.01),]
#df.plot.p$enrichment = rev(-log10(as.numeric(df.plot.p$elimKS)))
df.plot.p$enrichment = -log10(as.numeric(df.plot.p$elimKS))
colnames(df.plot.p)[2] <- "term"
colnames(df.plot.p)[7] <- "process"

#add separator bars
df.plot.p[14,] = c("","none1","","","","","none1",1)
df.plot.p[15,] = c("","none2","","","","","none2",1)
df.plot.p$enrichment = as.numeric(df.plot.p$enrichment)


df.plot.p$procnum = ifelse(df.plot.p$process == "bp",1,
                           ifelse(df.plot.p$process == "none1",2,
                                  ifelse(df.plot.p$process == "cc",3,
                                         ifelse(df.plot.p$process == "none2",4,5))))

b=df.plot.p %>%
  mutate(reord = as.numeric(procnum) + as.numeric(enrichment),
         term = fct_reorder(term, reord, .desc = F)) %>% 
  ggplot(aes(y =  reorder(term, -procnum),  x = enrichment, fill = process)) +
  geom_bar(stat = "identity",  width = 0.8,alpha=0.75)+
  scale_fill_manual(values = c("bp" = "#D59935","cc" = "#2A5B6E", "mf" = "#BE3B38"),na.value = NA)+
  theme_minimal()+
  ylab("")+
  xlab("")+
  ggtitle("PAR summer")+
  theme(legend.position = "none")+
  xlim(0,3.3)

#kd490sp
load("kd490spring_sig.gos.cats.Rdata")

df.plot.k = kd490spring_sig.gos[which(kd490spring_sig.gos$elimKS <= 0.01),]
# df.plot.k$enrichment = rev(-log10(as.numeric(df.plot.k$elimKS)))
df.plot.k$enrichment = -log10(as.numeric(df.plot.k$elimKS))

colnames(df.plot.k)[2] <- "term"
colnames(df.plot.k)[7] <- "process"

#add separator bars
df.plot.k[7,] = c("","","","","","","none1",1)
#df.plot.k[15,] = c("","none2","","","","","none2",1)
df.plot.k$enrichment = as.numeric(df.plot.k$enrichment)


df.plot.k$procnum = ifelse(df.plot.k$process == "bp",1,
                           ifelse(df.plot.k$process == "none1",2,
                                  ifelse(df.plot.k$process == "cc",3,
                                         ifelse(df.plot.k$process == "none2",4,5))))


c=df.plot.k %>%
  mutate(reord = as.numeric(procnum) + as.numeric(enrichment),
         term = fct_reorder(term, reord, .desc = F)) %>% 
  ggplot(aes(y =  reorder(term, -procnum),  x = enrichment, fill = process)) +
  geom_bar(stat = "identity",  width = 0.8,alpha=0.75)+
  scale_fill_manual(values = c("bp" = "#D59935","mf" = "#BE3B38"),na.value = NA)+
  theme_minimal()+
  ylab("")+
  xlab("")+
  ggtitle("kd490 spring")+
  theme(legend.position = "none")+
  xlim(0,3.3)

library(ggplot2)
library(gridExtra)
library(cowplot)

plot_grid(a, b, c, align = "v", nrow = 3, rel_heights = c(7.13/28,12.75/28,8.12/28))


