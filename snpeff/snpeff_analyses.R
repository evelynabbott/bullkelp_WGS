
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(ggpubr)

rm(list=ls())

load("env_65.Rdata") #environmental data
load("impact.prot.Rdata") #variants matched with protein IDs and genome location
high = impact.prot[impact.prot$impact=="HIGH",-c(2)] %>% distinct()

highfreq = high %>% count(pos, region,name = 'count')
highfreq$pos = as.numeric(highfreq$pos)
highfreq = highfreq %>% distinct()


bkgd = data.frame("scaf" = seq(1,38),rep=27,"bkc" = rep(c("black","white"),19))

#HIGH
a=ggplot(data=NULL)+
  geom_col(data=bkgd, aes(scaf,rep,fill=bkc),alpha=0.3,width = 1,just=0)+
  scale_fill_manual(values=c("lightgray","darkgray"))+
  geom_point(data=highfreq,aes(pos,count,color=region),size=3,alpha=0.8)+
  geom_smooth(data=highfreq,aes(pos,count,group=region,color=region,))+
  scale_color_manual(values=c("darkorange3","hotpink2","orchid4"))+
  theme_classic()+
  #scale_x_continuous(limits = c(0.3,38.5), expand = c(0, 0),breaks=seq(0,38,1))+
  scale_x_continuous(limits = c(0.8,38), expand = c(0, 0),breaks=seq(0,38,1))+
  scale_y_continuous(limits = c(0,27), expand = c(0, 0),breaks=seq(0,38,5))+
  theme(axis.ticks = element_blank())+
  xlab("scaffold")+
  ylab("variants")+
  theme(legend.position = "none",
        panel.background = element_rect(fill='transparent'))


Env.dataP = env #make names match
high = impact.prot[impact.prot$impact=="HIGH",-c(2)] %>% distinct()
high = high[,c(1,4)]
high = high[-which(high$pop =="FW"),]

popnum = data.frame(table(substr(rownames(Env.dataP),1,2)))
popnum = popnum[-which(popnum$Var1 == "FW"),]
colnames(popnum)=c("pop","numind")
popnum$pop = as.character(popnum$pop)

highfreq1 = high %>% count(pos, pop,name = 'count')
freqpopnum = full_join(highfreq1,popnum,by="pop")

fixed = freqpopnum[which(freqpopnum$numind == freqpopnum$count),]

table(fixed$pop)

fixtotal = data.frame(table(fixed$pop))
str(fixtotal)
colnames(fixtotal) = c("pop","total")
fixtotal$pop = as.character(fixtotal$pop)
fixtotal = transform(fixtotal,variables = reorder(pop,-total))

Group1<-c("FW","MC","JT","AH","PT","KR","FB")
Group2<-c("DB","SH","PP","HI","ED","SM","EB")
Group3<-c("LP","PV","SB","DI","SQ")

GenGroupK3<-1:nrow(fixtotal)
GenGroupK3[fixtotal$pop%in%Group1]<-"SJF"
GenGroupK3[fixtotal$pop%in%Group2]<-"WB"
GenGroupK3[fixtotal$pop%in%Group3]<-"SPS"

fixtotal$region = GenGroupK3

b=ggplot(fixtotal,aes(variables,total,fill=region))+
  geom_bar(stat="identity",alpha=0.85)+
  scale_fill_manual(values = c("hotpink2","darkorange3","orchid4"))+
  theme_minimal()+
  xlab("")+
  ylab("fixed variants")+
  theme(legend.position = "none")


grid.arrange(a,b,nrow=2)





#who has these variants? ----------------
#where are they?
rm(list=ls())
load("impact.prot.Rdata")

no.prot = impact.prot[,-c(2)] %>%  add_column(chr = gsub("\\..*","",impact.prot$pos)) %>% distinct() 
h.no.prot = no.prot[no.prot$impact == "HIGH",]
length(unique(h.no.prot$chr))
length(unique(h.no.prot$pos))


bypos = h.no.prot %>% group_split(pos)
byregion = h.no.prot %>% group_split(region)
byind = h.no.prot %>% group_split(ind)
bypop = h.no.prot %>% group_split(pop)

sapply(bypop, nrow)



posls = vector(length=37)
for (i in seq_along(bypos)) {
  posls[i] = unique(bypos[[i]][1])
}

posls = unlist(posls)
posls = unique(posls)

indperpos = data.frame("pos" = posls,"indnum" = sapply(bypos, nrow))
mean(indperpos$indnum)
ggplot(indperpos,aes(sample=indnum))+
  stat_qq()

#unique loci to population?
rls = list()
for (i in seq_along(bypos)) {
  rls[[i]] = unique(bypos[[i]][5])
}

byreg = data.frame("reg"=sapply(rls, nrow))


shapiro.test(indperpos$indnum)
sd(indperpos$indnum)

tab4supp = data.frame("variant_position" = indperpos$pos,"individuals" = indperpos$indnum, "regions" = byreg$reg)
#write.table(tab4supp,file = "variants_shared.csv",sep = ",",quote = F,col.names = T,row.names = F)

#annotations ---------
# length(unique(impact.prot.high$proteinId))
# GO.DF<-read.delim("Nerluet1_FilteredModels1_go_2024-03-15.tab")
# 
# GO.DF.snpeff = GO.DF[which(GO.DF$proteinId %in% impact.prot.high$proteinId),]
# length(unique(impact.prot.high$proteinId))
# length(unique(impact.prot.high$pos))
# length(unique(wb$proteinId))
# 
# Nereo.gff<-read.delim("Nerluet1_FilteredModels1_2024-03-15.gff3",header=F,skip=2)
# Nereo.gff$proteinId<-as.numeric(str_extract(Nereo.gff$V9,regex("(?<=(proteinId=))[0-9]+")))
# Nereo.gff = na.omit(Nereo.gff)
# 
# length(which(unique(impact.prot.high$proteinId) %in% Nereo.gff$proteinId))
# length(which(Nereo.gff$proteinId %in% GO.DF$proteinId))
# 
# GO.DF.snpeff$go_name
# [1] "nucleic acid binding"                                                                           
# [2] "helicase activity"                                                                              
# [3] "ATP binding"                                                                                    
# [4] "nucleic acid binding"                                                                           
# [5] "methyltransferase activity"                                                                     
# [6] "methylation"                                                                                    
# [7] "alkaline phosphatase activity"                                                                  
# [8] "glucuronosyl-N-acetylglucosaminyl-proteoglycan 4-alpha-N-acetylglucosaminyltransferase activity"
# [9] "N-acetylglucosaminyl-proteoglycan 4-beta-glucuronosyltransferase activity"                      
# [10] "nucleic acid binding"



#new box plots -------------------------------------------
rm(list=ls())
load("impact.prot.Rdata")
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(gridExtra)

impact.prot$region = factor(impact.prot$region,levels=c("PS","WI","SJF"))
impact.prot$ind = as.character(impact.prot$ind)

impact.prot = impact.prot[,-c(2)] %>% distinct()


grp = impact.prot[impact.prot$impact == "HIGH",]
grpnum = as.data.frame(table(grp$ind))
colnames(grpnum) = c("ind","num")
grpnum$ind = as.character(grpnum$ind)
grpregion = grp[,c("ind","region")] %>% distinct()
high = full_join(grpnum,grpregion,by="ind") %>% add_column(var = "high")

grp = impact.prot[impact.prot$impact == "MODERATE",]
grpnum = as.data.frame(table(grp$ind))
colnames(grpnum) = c("ind","num")
grpnum$ind = as.character(grpnum$ind)
grpregion = grp[,c("ind","region")] %>% distinct()
mod = full_join(grpnum,grpregion,by="ind") %>% add_column(var = "mod")

grp = impact.prot[impact.prot$impact == "LOW",]
grpnum = as.data.frame(table(grp$ind))
colnames(grpnum) = c("ind","num")
grpnum$ind = as.character(grpnum$ind)
grpregion = grp[,c("ind","region")] %>% distinct()
low = full_join(grpnum,grpregion,by="ind") %>% add_column(var = "low")

#normal distribution?
ggdensity(high$num)
ggqqplot(high$num)

shapiro.test(high$num) #normal
shapiro.test(mod$num) #not normal
shapiro.test(low$num) #not normal

high.stat=compare_means(num ~ region,  data = high,method = "anova")
mod.stat=compare_means(num ~ region,  data = mod)


my_comparisons <- list( c("PS", "SJF"), c("PS", "WI"), c("WI", "SJF") )

a=ggboxplot(high, x = "region", y = "num",fill="region",alpha=0.7)+
  stat_compare_means(comparisons = my_comparisons)+
  stat_compare_means(label.y = 26,label.x =1.1)+
  scale_fill_manual(values = rev(c("hotpink2","orchid4","darkorange3")))+
  ggtitle("high")+
  theme(legend.position = "none")+
  ylab("")+
  scale_x_discrete(labels=c("SPS","WB","SJF"))+
  xlab("")

b=ggboxplot(mod, x = "region", y = "num",fill="region",alpha=0.7)+
  stat_compare_means(comparisons = my_comparisons)+
  stat_compare_means(label.y=450, label.x = 1.1)+
  scale_fill_manual(values = rev(c("hotpink2","orchid4","darkorange3")))+
  scale_y_continuous(breaks = function(x) unique(floor(pretty(seq(min(x), (max(x) + 1) * 1.1)))))+
  ggtitle("moderate")+
  theme(legend.position = "none")+
  ylab("")+
  scale_x_discrete(labels=c("SPS","WB","SJF"))+
  xlab("")



c=ggboxplot(low, x = "region", y = "num",fill="region",alpha=0.7)+
  stat_compare_means(comparisons = my_comparisons)+
  stat_compare_means(label.y=340 ,label.x = 1.1)+
  scale_fill_manual(values = rev(c("hotpink2","orchid4","darkorange3")))+
  scale_y_continuous(breaks = function(x) unique(floor(pretty(seq(min(x), (max(x) + 1) * 1.1)))))+
  ggtitle("low")+
  theme(legend.position = "none")+
  ylab("")+
  scale_x_discrete(labels=c("SPS","WB","SJF"))+
  xlab("")


grid.arrange(a,b,c,nrow=1)


#stats
#high
mean(high[high$region =="PS",]$num)  #13.94444
mean(high[high$region =="WI",]$num)  #9.964286
mean(high[high$region =="SJF",]$num) #11.29412

max(high[high$region =="PS",]$num)  #19
max(high[high$region =="WI",]$num)  #18
max(high[high$region =="SJF",]$num) #20

min(high[high$region =="PS",]$num)  #7
min(high[high$region =="WI",]$num)  #5
min(high[high$region =="SJF",]$num) #5

#mod
mean(mod[mod$region =="PS",]$num)  #330.7778
mean(mod[mod$region =="WI",]$num)  #341.3571
mean(mod[mod$region =="SJF",]$num) #306.8235

max(mod[mod$region =="PS",]$num)  #384
max(mod[mod$region =="WI",]$num)  #398
max(mod[mod$region =="SJF",]$num) #353

min(mod[mod$region =="PS",]$num)  #304
min(mod[mod$region =="WI",]$num)  #307
min(mod[mod$region =="SJF",]$num) #235


#low
mean(low[low$region =="PS",]$num)  #263.5556
mean(low[low$region =="WI",]$num)  #267.8571
mean(low[low$region =="SJF",]$num) #249.0588

max(low[low$region =="PS",]$num)  #288
max(low[low$region =="WI",]$num)  #304
max(low[low$region =="SJF",]$num) #283

min(low[low$region =="PS",]$num)  #236
min(low[low$region =="WI",]$num)  #244
min(low[low$region =="SJF",]$num) #183

min(low$num)

hist(low$num,breaks=350)



################################# --------------------------------------------
##### MISSENSE SYNONYMOUS #######
#################################

rm(list=ls())
load("envdatap.Rdata")

#box plot first
load("snpeff_ind_miss_syn.Rdata")
my_comparisons <- list( c("SPS", "SJF"), c("SPS", "WB"), c("WB", "SJF") )
missvars_by_ind$region = factor(missvars_by_ind$region,levels=c("SPS","WB","SJF"))
synvars_by_ind$region = factor(synvars_by_ind$region,levels=c("SPS","WB","SJF"))





ggboxplot(missvars_by_ind,x="region",y="var")+
  stat_compare_means(comparisons = my_comparisons)

ggboxplot(synvars_by_ind,x="region",y="var")+
  stat_compare_means(comparisons = my_comparisons)

#curves
rm(list=ls())
load("snpeff_miss_syn.Rdata")
#some exploration
barplot(prop.table(table(synvars$site)))

popnum = data.frame(table(substr(rownames(Env.dataP),1,2)))
colnames(popnum)=c("pop","freq")

popcount = synvars %>% group_by(pos,site)%>% summarize ("Freq" = n())
popcount = full_join(popnum,)


# table(missvars$pos)
# 
# synvars$indnum = as.numeric(factor(synvars$ind))
# 
# syn.plot=synvars %>%
#   count(pos, region,name = 'count')
# 
# miss.plot=missvars %>%
#   count(pos, region,name = 'count')
# 
# 
# ggplot(syn.plot,aes(pos,count,group = region,color=region)) +
#   geom_smooth(linewidth = 2)+
#   scale_color_manual(values = c("hotpink2","orchid4","darkorange3"))
#   
# ggplot(miss.plot,aes(pos,count,group = region,color=region)) +
#   geom_smooth(linewidth = 2)+
#   scale_color_manual(values = c("hotpink2","orchid4","darkorange3"))

sjf.syn = synvars[synvars$region == "SJF",]
wb.syn = synvars[synvars$region == "WI",]
sps.syn = synvars[synvars$region == "PS",]

sjf.syn=sjf.syn %>% group_by(pos)%>% summarize ("Freq" = n()) %>% add_column(region="SJF")
wb.syn=wb.syn %>% group_by(pos)%>% summarize ("Freq" = n()) %>% add_column(region="WB")
sps.syn=sps.syn %>% group_by(pos)%>% summarize ("Freq" = n()) %>% add_column(region="SPS")

any = synvars %>% count(pos,name="count") 

highfreq1=rbind(sjf.syn,wb.syn,sps.syn)
all$prop = all$Freq/63

ggplot(all,aes(pos,Freq,group = region,color=region)) +
  geom_smooth(linewidth = 2)+
  geom_point(alpha=0.2)+
  geom_smooth



