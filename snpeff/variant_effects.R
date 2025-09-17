
#make tables showing types of mutations for each variant
rm(list=ls())

files <- list.files(pattern = "*ANN")
load("snpeff_results.Rdata")
test=filesepd[[1]]

subsetcols = lapply(filesepd, function(x) subset(x, select = c(1:7)))

subsetcols <- subsetcols %>%
  map(~ .x %>%
        unite(pos, `column_a`,
              `column_b`,remove = TRUE, sep="."))

#add names
subsetcols = mapply(cbind, subsetcols, "ind"=gsub("\\_.*","",files), SIMPLIFY=F)

all = bind_rows(subsetcols)

Group1<-c("FW","MC","JT","AH","KR")
Group2<-c("DB","SH","PP","HI","ED","SM","EB","PT","FB")
Group3<-c("LP","PV","SB","DI","SQ")

#The Vector GenGroupK3 will have the membership code for each site to one of the three genetic groups
all$pop = ifelse(substr(all$ind,1,2) %in% Group1, "SJF",
                 ifelse(substr(all$ind,1,2) %in% Group2, "WB", "SPS"))


#HIGH ----------------
htot=data.frame(table(all[all$column_g == "HIGH",]$column_f))
colnames(htot) = c("Effect","Total")

hsjf=data.frame(table(all[all$column_g == "HIGH" & all$pop == "SJF",]$column_f))
colnames(hsjf) = c("Effect","SJF")

hwb=data.frame(table(all[all$column_g == "HIGH" & all$pop == "WB",]$column_f))
colnames(hwb) = c("Effect","WB")

hsps=data.frame(table(all[all$column_g == "HIGH" & all$pop == "SPS",]$column_f))
colnames(hsps) = c("Effect","SPS")


htab = full_join(hsps, full_join(hwb, full_join(hsjf,htot)))
htab$Prop = htab$Total/sum(htab$Total) * 100

#MODERATE ----------------
mtot=data.frame(table(all[all$column_g == "MODERATE",]$column_f))
colnames(mtot) = c("Effect","Total")

msjf=data.frame(table(all[all$column_g == "MODERATE" & all$pop == "SJF",]$column_f))
colnames(msjf) = c("Effect","SJF")

mwb=data.frame(table(all[all$column_g == "MODERATE" & all$pop == "WB",]$column_f))
colnames(mwb) = c("Effect","WB")

msps=data.frame(table(all[all$column_g == "MODERATE" & all$pop == "SPS",]$column_f))
colnames(msps) = c("Effect","SPS")

mtab = full_join(msps, full_join(mwb, full_join(msjf,mtot)))
mtab$Prop = signif(mtab$Total/sum(mtab$Total),digits=2)


#LOW ----------------
ltot=data.frame(table(all[all$column_g == "LOW",]$column_f))
colnames(ltot) = c("Effect","Total")

lsjf=data.frame(table(all[all$column_g == "LOW" & all$pop == "SJF",]$column_f))
colnames(lsjf) = c("Effect","SJF")

lwb=data.frame(table(all[all$column_g == "LOW" & all$pop == "WB",]$column_f))
colnames(lwb) = c("Effect","WB")

lsps=data.frame(table(all[all$column_g == "LOW" & all$pop == "SPS",]$column_f))
colnames(lsps) = c("Effect","SPS")

ltab = full_join(lsps, full_join(lwb, full_join(lsjf,ltot)))
#ltab$Prop = signif(ltab$Total/sum(ltab$Total),digits=2)
ltab$Prop = ltab$Total/sum(ltab$Total) * 100


write.table(htab,file="snpeff_high_effect.csv",quote = F,col.names = T,row.names = F,sep = ",")
write.table(mtab,file="snpeff_moderate_effect.csv",quote = F,col.names = T,row.names = F,sep = ",")
write.table(ltab,file="snpeff_low_effect.csv",quote = F,col.names = T,row.names = F,sep = ",")


#compare effects to lfmm
all2comp = all
all2comp$pos = gsub(".*_","",all2comp$pos)

comps = all2comp[which(all2comp$pos %in% pos4lfmm),]
comps = comps[comps$column_g != "MODIFIER",]

table(comps$column_g)
small = comps[,c(1,5,6)] %>% distinct()





