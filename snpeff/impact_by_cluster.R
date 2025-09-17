
#SNPEff - initial exploration
library(tidyverse)
library(stringr)
library(gridExtra)
library(ggpubr)


rm(list=ls())
load("filesepd.Rdata")
files <- list.files(pattern = "*ANN") 

#look at impact level ------------------------------
impactlvl = lapply(filesepd, function(x) subset(x, select = c(1,2,7)))
testdf=impactlvl[[1]]

lowimp = lapply(impactlvl, function(x) subset(x, column_g == "LOW"))
modimp = lapply(impactlvl, function(x) subset(x, column_g == "MODERATE"))
highimp = lapply(impactlvl, function(x) subset(x, column_g == "HIGH"))
modifimp = lapply(impactlvl, function(x) subset(x, column_g == "MODIFIER"))

#save(list = ls(pattern = '.imp'), file = "snp_impact.RData")

low_by_ind=sapply(lowimp, nrow)
mod_by_ind=sapply(modimp, nrow)
high_by_ind=sapply(highimp, nrow)
modif_by_ind=sapply(modifimp, nrow)

ind_names = files

pops = data.frame("pop"=substr(ind_names,1,2))

Group1<-c("FW","MC","JT","AH","PT","KR","FB")
Group2<-c("DB","SH","PP","HI","ED","SM","EB")
Group3<-c("LP","PV","SB","DI","SQ")

GenGroupK3<-1:length(ind_names)
GenGroupK3[pops$pop%in%Group1]<-"SJF"
GenGroupK3[pops$pop%in%Group2]<-"WI"
GenGroupK3[pops$pop%in%Group3]<-"PS"

pops$region = GenGroupK3

low_by_ind=data.frame("vars"=sapply(lowimp, nrow),"ind"=ind_names,"pop" = pops$pop,"region" = pops$region)
mod_by_ind=data.frame("vars"=sapply(modimp, nrow),"ind"=ind_names,"pop" = pops$pop,"region" = pops$region)
high_by_ind=data.frame("vars"=sapply(highimp, nrow),"ind"=ind_names,"pop" = pops$pop,"region" = pops$region)
modif_by_ind=data.frame("vars"=sapply(modifimp, nrow),"ind"=ind_names,"pop" = pops$pop,"region" = pops$region)


my_comparisons <- list( c("PS", "SJF"), c("PS", "WI"), c("WI", "SJF") )

newlab = c("SPS","WB","SJF")
low_by_ind$region = factor(low_by_ind$region, levels = c("PS","WI","SJF"))
mod_by_ind$region = factor(mod_by_ind$region, levels = c("PS","WI","SJF"))
high_by_ind$region = factor(high_by_ind$region, levels = c("PS","WI","SJF"))


a=ggboxplot(low_by_ind, x = "region", y = "vars",fill="region",alpha=0.7)+
  stat_compare_means(comparisons = my_comparisons)+
  scale_fill_manual(values = rev(c("hotpink2","orchid4","darkorange3")))+
  theme(legend.position="none")+
  ggtitle("low")+
  scale_x_discrete(labels= newlab)+
  ylab("number of variants")

b=ggboxplot(mod_by_ind, x = "region", y = "vars",fill="region",alpha=0.7)+
  stat_compare_means(comparisons = my_comparisons)+
  scale_fill_manual(values = rev(c("hotpink2","orchid4","darkorange3")))+
  theme(legend.position="none")+
  ggtitle("moderate")+
  scale_x_discrete(labels= newlab)+
  ylab("number of variants")

c=ggboxplot(high_by_ind, x = "region", y = "vars",fill="region",alpha=0.7)+
  stat_compare_means(comparisons = my_comparisons)+
  scale_fill_manual(values = rev(c("hotpink2","orchid4","darkorange3")))+
  theme(legend.position="none")+
  ggtitle("high")+
  scale_x_discrete(labels= newlab)+
  ylab("number of variants")

grid.arrange(a,b,c,nrow=1)


#stats for manuscript
lowregion = low_by_ind[low_by_ind$region == "WI",]
mean(lowregion$vars)
min(lowregion$vars)
max(lowregion$vars)

impactlvl[1]
#make a plot: x = position, y = ind colored by l/m/h
#3 plots, one for each region. stacked.

test.df = as.data.frame(impactlvl[1])
#code position as scaff #.pos. SCAF_1-100 = 1.100

test.df = data.frame("pos"=paste0(sub('SCAF_',"",test.df$column_a),".",test.df$column_b),
                     "ind"=c(rep(1,times=nrow(test.df))),
                     "impact" = test.df$column_g)

test.df2 = as.data.frame(impactlvl[5])
test.df2 = data.frame("pos"=paste0(sub('SCAF_',"",test.df2$column_a),".",test.df2$column_b),
                      "ind"=c(rep(2,times=nrow(test.df2))),
                      "impact" = test.df2$column_g)

test.plot = rbind(test.df,test.df2)

#get rid of modifier
#`%!in%` = Negate(`%in%`)

test.nomodif = test.plot[which(test.plot$impact != "MODIFIER"),]

ggplot(test.nomodif,aes(pos,ind,color=impact))+
  geom_point(alpha=0.2)+
  scale_color_manual(values=c("red","grey","grey"))+
  ylim(1,2)

#try on entire list -----------------------
impactlvl.list = impactlvl
ID = seq(1,63,by=1)

pops$ID = files

#assign sample order by region, quick n dirty - for the plot
poporder = pops
poporderps = poporder[which(poporder$region == "PS"),]
poporderps$order = seq(1,nrow(poporderps),by = 1)
poporderwi = poporder[which(poporder$region == "WI"),]
poporderwi$order = seq(19,46,by = 1)
popordersjf = poporder[which(poporder$region == "SJF"),]
popordersjf$order = seq(47,63,by = 1)
poporder = rbind(poporderps,poporderwi,popordersjf)
pops$originalorder = c(1:63)

pops = right_join(pops,poporder,by="ID")
pops = data.frame("pop" = pops$pop.x,"region" = pops$region.x,"order" = pops$order)


for( i in seq_along(impactlvl.list)){
  impactlvl.list[[i]] <- data.frame("ind" = rep(ID[i],nrow(impactlvl.list[[i]])),
                                    "pos" = paste0(sub('SCAF_',"",impactlvl.list[[i]]$column_a),".",impactlvl.list[[i]]$column_b),
                                    "impact" = impactlvl.list[[i]]$column_g,
                                    "pop" = pops$pop [i],
                                    "region" = pops$region [i],
                                    "order" = pops$order [i])
}

#rbind into one massive dataframe
df.all=bind_rows(impactlvl.list)

impact.all = df.all
#save(impact.all, file = "impact_snpbyind.Rdata")


