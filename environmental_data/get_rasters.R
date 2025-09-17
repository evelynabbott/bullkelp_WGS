
# sourcing environmental data from NOAA coastwatch database
# http://geog.uoregon.edu/bartlein/courses/geog490/week04-netCDF.html

library(raster)
#install.packages("ncdf4")
library(ncdf4)

#ABBREVIATION LEGEND 
# SST : seasurface temp
# POC : Particulate organic carbon
# PIC : Pollution Identification and Correction
# ChlA : Chlorophyl A
# PAR: coefficient for photosynthetically active radiation

#################################### ---------------------------------------
###### DOWNLOAD .nc FILES ##########
####################################


#kd490
# urlq="https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdMH1kd490mday_R2022SQ.nc?Kd_490%5B(2003-01-16):1:(2022-02-16T00:00:00Z)%5D%5B(47):1:(48.3)%5D%5B(-124.3):1:(-122.28)%5D"
# urlq= paste0("\"",urlq,"\"")
# system2(command="curl", args= c("-g",urlq,"-o","kd490.PugeSound.2012-2023.nc"))  #checked
# kd490<-stack("kd490.PugeSound.2012-2023.nc")    #This ends up being the one used

# 
# #SST
# urlq="https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdMBsstdmday.nc?sst%5B(2006-10-17T12:00:00Z):1:(2023-10-16T12:00:00Z)%5D%5B(0.0):1:(0.0)%5D%5B(47):1:(48.3)%5D%5B(235.6):1:(237.92)%5D"
# urlq= paste0("\"",urlq,"\"")
# system2(command="curl", args= c("-g",urlq,"-o","SST.PugeSound.2012-2023.nc"))  #checked
# SST<-stack("SST.PugeSound.2012-2023.nc")    #This ends up being the one used

# 
# 
# 
#SST monthly NOAA Global Coral Bleaching Monitoring, 5km NOAA Global Coral Bleaching Monitoring, 5km See license in site https://coastwatch.pfeg.noaa.gov/erddap/griddap/NOAA_DHW_monthly.html
#See also here for full details https://coralreefwatch.noaa.gov/product/5km/methodology.php#sst
# urlq="https://coastwatch.pfeg.noaa.gov/erddap/griddap/NOAA_DHW_monthly.nc?sea_surface_temperature%5B(1985-01-16):1:(2023-10-16T00:00:00Z)%5D%5B(47):1:(49)%5D%5B(-124.3):1:(-122.28)%5D,mask%5B(1985-01-16):1:(2023-10-16T00:00:00Z)%5D%5B(47):1:(49)%5D%5B(-124.3):1:(-122.28)%5D"
# urlq= paste0("\"",urlq,"\"")
# system2(command="curl", args= c("-g",urlq,"-o","SSTmonthly.PugeSound.2010-2020.nc")) #will need to add the southernmost sites from buoy data
# SST.month<-stack("SSTmonthly.PugeSound.2010-2020.nc")
# 
# # 
# # #PAR     https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdMWpar0mday.html
# # #PAR ( a second try which seems to have great cover https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdVH2018parmday.html
# # 
# #This is PAR2nd below
# urlq="https://coastwatch.pfeg.noaa.gov/erddap/griddap/erdVH2018parmday.nc?par%5B(2012-01-15):1:(2022-07-15T00:00:00Z)%5D%5B(47):1:(48.3)%5D%5B(-124.5):1:(-122)%5D"
# urlq= paste0("\"",urlq,"\"")
# system2(command="curl", args= c("-g",urlq,"-o","PAR2nd.PugeSound.2002-2023.nc"))  #checked
# PAR2nd<-stack("PAR2nd.PugeSound.2002-2023.nc")

# 
# 
# #kdPAR PAr attenuation, it might be interesting just recent years experimental but cover the puget SOunds
# #no work
# #monthly https://coastwatch.pfeg.noaa.gov/erddap/griddap/nesdisVHNSQkdparMonthly.html
# #urlq="https://coastwatch.pfeg.noaa.gov/erddap/griddap/nesdisVHNSQkdparMonthly.nc?kd_par%5B(2012-01-02T12:00:00Z):1:(2023-10-02T12:00:00Z)%5D%5B(0.0):1:(0.0)%5D%5B(47):1:(48.3)%5D%5B(-124.5):1:(-122)%5D"
# #urlq= paste0("\"",urlq,"\"")
# 

#SalinityIPRC<-stack("SalinityAquar.PugeSound.2012-2022.nc")
 

#save(PAR2nd,PIC,POC,SST,SST.month,kd490,SalinityIPRC, file = "rasters.Rdata")


#save(PAR2nd,SST.month,kd490,SalinityIPRC, file = "rasters2impute.Rdata")


# load data ----
rm(list=ls())
#load("rasters.Rdata")

load("rasters2impute.Rdata")

#POC plots ------------------------------------------------------------

#Salinity -------------------------------
MonthsSAL<-as.numeric(substr(names(SalinityIPRC),7,9))

cuts=c(seq(266.6,2000,100), seq(2500,13000,2000))#set breaks
pal <- colorRampPalette(c("skyblue1","seagreen2","yellow","red4"))

#Seasonal averages, NOTE the seasonal rasters produced
Winter.Avg.salinity<-calc(SalinityIPRC[[which(MonthsSAL%in%c(12,1,2))]],function(x){mean(x, na.rm=T)})
plot(Winter.Avg.salinity,main="Winter SAL",breaks=cuts, col = pal(length(cuts)))

Spring.Avg.sal<-calc(SalinityIPRC[[which(MonthsSAL%in%c(12,1,2))]],function(x){mean(x, na.rm=T)})
plot(Spring.Avg.sal,main="Spring POC",breaks=cuts, col = pal(length(cuts)))

Summmer.Avg.sal<-calc(SalinityIPRC[[which(MonthsSAL%in%c(12,1,2))]],function(x){mean(x, na.rm=T)})
plot(Summmer.Avg.sal,main="Summer POC",breaks=cuts, col = pal(length(cuts)))

Fall.Avg.sal<-calc(SalinityIPRC[[which(MonthsSAL%in%c(12,1,2))]],function(x){mean(x, na.rm=T)})
plot(Fall.Avg.sal,main="Fall POC",breaks=cuts, col = pal(length(cuts)))

#save(Winter.Avg.salinity,Spring.Avg.sal,Summmer.Avg.sal,Fall.Avg.sal,file = "Salinity_rasters.Rdata")


#2ndPAR ------------------------------------------
MonthsPAR<-as.numeric(substr(names(PAR2nd),7,9))

cuts=c(seq(7,51,1))#set breaks


#Seasonal averages, NOTE the seasonal rasters produced
Winter.Avg.PAR2<-calc(PAR2nd[[which(MonthsPAR%in%c(1,2))]],function(x){mean(x, na.rm=T)})
plot(Winter.Avg.PAR2,main="Winter PAR",breaks=cuts,col = pal(length(cuts)))

Spring.Avg.PAR2<-calc(PAR2nd[[which(MonthsPAR%in%c(3,4,5))]],function(x){mean(x, na.rm=T)})
plot(Spring.Avg.PAR2,main="Spring PAR",breaks=cuts,col = pal(length(cuts)))

Summmer.Avg.PAR2<-calc(PAR2nd[[which(MonthsPAR%in%c(6,7,8))]],function(x){mean(x, na.rm=T)})
plot(Summmer.Avg.PAR2,main="Summer PAR",breaks=cuts,col = pal(length(cuts)))

Fall.Avg.PAR2<-calc(PAR2nd[[which(MonthsPAR%in%c(9,10,11))]],function(x){mean(x, na.rm=T)})
plot(Fall.Avg.PAR2,main="Fall PAR",breaks=cuts,col = pal(length(cuts)))

#save(PAR2nd,Winter.Avg.PAR2,Spring.Avg.PAR2,Summmer.Avg.PAR2,Fall.Avg.PAR2,file="PAR_rasters.Rdata")

#kd490 ----------------------------------
#Seasonal averages, NOTE the seasonal rasters produced
plot(kd490)
MonthsPAR<-as.numeric(substr(names(kd490),7,9))
#cuts=c(seq(7,51,1))#set breaks
pal = colorRampPalette(rev(c("#fde725", "#5ec962","#21918c", "#3b528b","#440154")))

Winter.Avg.kd490<-calc(kd490[[which(MonthsPAR%in%c(1,2))]],function(x){mean(x, na.rm=T)})
plot(Winter.Avg.kd490,main="Winter PAR",col = pal(length(cuts)))

Spring.Avg.kd490<-calc(kd490[[which(MonthsPAR%in%c(3,4,5))]],function(x){mean(x, na.rm=T)})
plot(Spring.Avg.kd490,main="Spring PAR",col = pal(length(cuts)))

Summmer.Avg.kd490<-calc(kd490[[which(MonthsPAR%in%c(6,7,8))]],function(x){mean(x, na.rm=T)})
plot(Summmer.Avg.kd490,main="Summer PAR",col = pal(length(cuts)))

Fall.Avg.kd490<-calc(kd490[[which(MonthsPAR%in%c(9,10,11))]],function(x){mean(x, na.rm=T)})
plot(Fall.Avg.kd490,main="Fall PAR",col = pal(length(cuts)))

#save(kd490,Winter.Avg.kd490,Spring.Avg.kd490,Summmer.Avg.kd490,Fall.Avg.kd490,file="kd490_rasters.Rdata")


# SST PLOTS -----------------------------
names(SST.month)
MonthsSSTmonth<-as.numeric(substr(names(SST.month),7,9))

cuts=c(seq(7,18,1))#set breaks
pal <- colorRampPalette(c("skyblue1","seagreen2","red4"))

par(mfrow=c(2,2))
Winter.Avg.SST.month<-calc(SST.month[[which(MonthsSSTmonth%in%c(12,1,2))]],function(x){mean(x, na.rm=T)})
plot(Winter.Avg.SST.month,main="Winter SST", breaks=cuts, col = pal(length(cuts)))

Spring.Avg.SST.month<-calc(SST.month[[which(MonthsSSTmonth%in%c(3,4,5))]],function(x){mean(x, na.rm=T)})
plot(Spring.Avg.SST.month,main="Spring SST",col = pal(length(cuts)))

Summmer.Avg.SST.month<-calc(SST.month[[which(MonthsSSTmonth%in%c(6,7,8))]],mean, na.rm=TRUE)
plot(Summmer.Avg.SST.month,main="Summer SST", col = pal(length(cuts)))


Fall.Avg.SST.month<-calc(SST.month[[which(MonthsSSTmonth%in%c(9,10,11))]],function(x){mean(x, na.rm=T)})
plot(Fall.Avg.SST.month,main="Fall SST",col = pal(length(cuts)))


#save(Winter.Avg.SST.month,Spring.Avg.SST.month,Summmer.Avg.SST.month,Fall.Avg.SST.month, file="SST_rasters.Rdata")

