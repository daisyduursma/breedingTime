rm(list = ls())
library(raster)
library(circular)
library(gtools)
library(Hmisc)
library(plyr)
library(Hmisc)
library(multcomp)
library(reshape2)
#figure
library(lme4)
#library(afex)
library(ENmisc)
library(car)
library(maptools)
library(doBy)


#working directory
dat.dir<-'/Users/daisy/Google Drive/PhD/Data/Observaitons/Cleaned/Breeding/'
#read in observations
dat<-read.csv(paste0(dat.dir,'PointOfLayDayOfYear2016-09-20.csv'))
#add koeppen zones
koeppen<-raster(paste0('/Users/daisy/Google Drive/PhD/Data/Spatial/BOM_climate_zones/kpngrp_major/koepenReclassified.asc'),
                proj4string=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))

aus<-readShapeSpatial("/Users/daisy/Google Drive/PhD/Data/Spatial/AustraliaPolygon/STE11aAust.shp",
                      proj4string=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))
CapitalTerr <- rasterize(aus[aus@data$STATE_CODE == 8, ],koeppen)
dat$ACT<-extract(CapitalTerr,data.frame(cbind(dat$lon, dat$lat)))
ACT<-subset(dat,ACT==1)
#convert data to radians to make circular vector
ACT$Radians <-(ACT$DOY_PL/366*360)*pi / 180
#traits for species
inc<-read.csv('/Users/daisy/Google Drive/PhD/BreedingTiming/tables/SpeciesOfInterest_2016-10-07.csv')
ACT$month<-formatC(ACT$month,width=2,flag='0') #format months 1 becomes 01
ACT$date<-paste0(ACT$year,ACT$month) #make date string that matches NOAA SOI
ACT1966<-subset(ACT, date >=196601) #limit time to same as Crick data
ACT1990<-subset(ACT, date >=199001) #limit time to same as Crick data

#get species with more than 100 Egg or young or multi observations from later than 1990
datAcc<-subset(ACT1990, type == "solitary-egg" | type == "multi"| type == "young")
datAcc_sp<-as.data.frame(table(datAcc$Scientific.Name))
datAcc_sp<-as.vector(droplevels(subset(datAcc_sp,Freq>=100)$Var1))
datAcc_sp<-data.frame(Species =datAcc_sp, group ="HighAcc")
datAcc<-datAcc[datAcc$Scientific.Name %in% datAcc_sp$Species,]
#get species with more than 100  observations from later than 1990
datACT2<-ACT1990[ACT1990$Scientific.Name %nin% datAcc_sp$Species,]
datACT2_sp<-as.data.frame(table(datACT2$Scientific.Name))
datACT2_sp<-as.vector(droplevels(subset(datACT2_sp,Freq>=100)$Var1))
datACT2<-datACT2[datACT2$Scientific.Name %in% datACT2_sp,]
fin_sp1990<-c(datACT2_sp,as.vector(datAcc_sp$Species))#list of species
ACT1990<-rbind(datAcc,datACT2)#data for later than 1990

#Data from 1966 to present
ACT1966<-ACT1966[ACT1966$Scientific.Name %nin% fin_sp1990,]
#get species with more than 100 Egg or young or multi observations from later than 1990
datAcc1966<-subset(ACT1966, type == "solitary-egg" | type == "multi"| type == "young")
datAcc1966_sp<-as.data.frame(table(datAcc1966$Scientific.Name))
datAcc1966_sp<-as.vector(droplevels(subset(datAcc1966_sp,Freq>=100)$Var1))
datAcc1966_sp<-data.frame(Species =datAcc1966_sp, group ="HighAcc")
datAcc1966<-datAcc1966[datAcc1966$Scientific.Name %in% datAcc1966_sp$Species,]
#get species with more than 50  observations from later than 1966
datACT1966_2<-ACT1966[ACT1966$Scientific.Name %nin% datAcc1966_sp$Species,]
datACT1966_2_sp<-as.data.frame(table(datACT1966_2$Scientific.Name))
datACT1966_2_sp<-as.vector(droplevels(subset(datACT1966_2_sp,Freq>=50)$Var1))
datACT1966_2<-datACT1966_2[datACT1966_2$Scientific.Name %in% datACT1966_2_sp,]


#bring the data back together
ACT<-droplevels(rbind(ACT1990,datAcc1966,datACT1966_2))
#list of species 
fin_sp<-as.vector(unique(ACT$Scientific.Name))

#make histograms of first egg dates
hist(ACT$DOY_PL,main="",xlab="",ylab = "",xlim = c(1,365),axes=FALSE,freq=FALSE,breaks=seq(1,366,length = 13))

ACTSummary<-list()
for (i in 1:length(fin_sp)){
    #get species data in specific phase
    spdat<-subset(ACT,Scientific.Name==fin_sp[i])#get species data
    NoObs<-as.numeric(nrow(spdat))#get observation count
    # CommonName<-as.character(subset(inc,ScientificName==fin_sp[i],"Common.me",drop=TRUE))
     Order<-as.character(subset(inc,ScientificName==fin_sp[i],order,drop=TRUE))
    # family<-as.character(subset(inc,ScientificName==fin_sp[i],family,drop=TRUE))
    # #clades<-as.character(subset(inc,ScientificName==fin_sp[i],"Patch.clade",drop=TRUE))
    obsRadians<-spdat$Radians
    #date of quantiles
    quant<-quantile.circular(obsRadians,c(.05,.5,.95),type=8)
    Per5<-round((quant[[1]]*180/pi)/360*365)
    Per50<-round((quant[[2]]*180/pi)/360*365)
    Per95<-round((quant[[3]]*180/pi)/360*365)
    #Breeding Period Length
    StartDate<-as.numeric(Per5)
    EndDate<-as.numeric(Per95)
    BPL<-ifelse (EndDate < StartDate, 365-StartDate+EndDate,EndDate - StartDate )
     # #put the data together and name it
    Perdat<-cbind(paste(fin_sp[i]),Order, NoObs,
                  Per5,Per50,Per95,BPL)
    colnames(Perdat)<-c('Species','Order','Observations','Quantile5',
                        'Quantile50','Quantile95','BreedingPeriod')
    ACTSummary[[i]]<-Perdat
}
  

ACTSummary<-as.data.frame(do.call('rbind', ACTSummary))


ACTSummary$BreedingPeriod <-as.numeric(levels(ACTSummary$BreedingPeriod))[ACTSummary$BreedingPeriod]
ACTSummary$Quantile5 <-as.numeric(levels(ACTSummary$Quantile5))[ACTSummary$Quantile5]
ACTSummary$Quantile50 <-as.numeric(levels(ACTSummary$Quantile50))[ACTSummary$Quantile50]
ACTSummary$Quantile95 <-as.numeric(levels(ACTSummary$Quantile95))[ACTSummary$Quantile95]



# Simple function for placing labels on a figure.
plotlabel <- function(txt, where, inset=0.1, font=2, inset.x=inset, inset.y=inset,...){
  u <- par()$usr
  if(grepl("left",where))x <- u[1] + inset.x*(u[2]-u[1])
  if(grepl("right",where))x <- u[2] - inset.x*(u[2]-u[1])
  if(grepl("bottom",where))y <- u[3] + inset.y*(u[4]-u[3])
  if(grepl("top",where))y <- u[4] - inset.y*(u[4]-u[3])
  
  text(x,y,txt,font=font,...)
}


#################
# What are are the differences in the NHTR vs. Temperate Australia? 
#################


ACTSummary<-subset(ACTSummary, select=c('Order','BreedingPeriod'))
ACTSummary$Region<-"Australia"

#how many had breeding periods more than 100 days long?
nrow(subset(ACTSummary,BreedingPeriod>=100))

#TEMPdat<- subset(BP, Region=="Temperate",select=c("BreedingPeriod","Region")) #AUS temp. breeding period length
UK<-subset(read.csv('/Users/daisy/Google Drive/PhD/BreedingTiming/tables/Joys_etal_Breeding_periods_England.csv')
           ,select=c('Order',"length.of.breeding")) #UK breeding period length
colnames(UK)<-c('Order',"BreedingPeriod")
UK$Region<-"UK"
df<-droplevels(rbind(ACTSummary,UK))#dataframe of UK and AUS temp

#Run ANOVA - mixed effect order as random

bp<-lmer(BreedingPeriod ~ Region + (1|Order),data=df)#add order
Anova(bp, test.statistic="F")


#means and SD of regions
summaryBy(BreedingPeriod ~ Region, data=df,
          FUN=c(mean,sd))



#################
# #transparent histograms of breeding period length, 1 column wide
#################
pdf(file = "/Users/daisy/Google Drive/PhD/BreedingTiming/figures/PeakBreedinPeriodUKandACT20161115.pdf",
    height = 4.5,
    width = 3.5)

par(mfrow=c(2,1), oma=c(4,4,5,0), mar=c(1,2,1,1), las=1, cex=.7, xaxs="i",
    yaxs = "i")
#set the intervals for the histograms
breaks <- seq(0,350,25)
#PREP DATA
d <- split(df, df$Region)
a = d[[1]][,1]
b = d[[2]][,1]
#UK
hist(b,freq=F,
     xlim=c(0,300), 
     col=rgb(0, 0, 0,0.5),
     breaks = breaks,
     main=NULL, 
     axes=FALSE,
     ylim=c(0,0.016))
plotlabel(expression(paste("NHTR (", bold("A"), ")")), "topright")

Axis(side=1, labels=FALSE,at=seq(0,350,25))
Axis(side=2, at=c(0,0.004,0.008,0.012,0.016),labels=c(0,0.004,0.008,0.012,0.016))
#AUS temperate
hist(a,freq=F,
     xlim=c(0,300), 
     col=rgb(0, 0, 0,0.5),
     breaks = breaks,
     main=NULL,
     axes=FALSE,
     ylim=c(0,0.016))
plotlabel(expression(paste("SHTR (", bold("B"), ")")), "topright")
Axis(side=1, labels=FALSE,at=seq(0,350,25))
Axis(side=1, labels=seq(0,350,50),at=seq(0,350,50))
Axis(side=2, at=c(0,0.004,0.008,0.012,0.016),labels=c(0,0.004,0.008,0.012,0.016))
#y and y labels
mtext("Frequency", 2, 2, outer=TRUE,las=0)
mtext("Breeding period (days)", 1, 2, outer=TRUE)

dev.off()





######################################
#what are average climatic conditions
#######################################

Eng<-readShapeSpatial('/Users/daisy/Google Drive/PhD/Data/Spatial/England/United_Kingdom_AL4.shp',
                      proj4string=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))

aus<-readShapeSpatial("/Users/daisy/Desktop/1259030001_ste11aaust_shape/STE11aAust.shp",
                      proj4string=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))
ACT<-aus[aus@data$STATE_CODE == 8, ]
CapitalTerr <- rasterize(aus[aus@data$STATE_CODE == 8, ],koeppen)

coldestQuart<-raster('/Users/daisy/Google Drive/PhD/Data/Spatial/Climate/WorldClim_bio_10m_bil/bio11.bil')
warmestQuart<-raster('/Users/daisy/Google Drive/PhD/Data/Spatial/Climate/WorldClim_bio_10m_bil/bio10.bil')


extract(coldestQuart,ACT,fun=mean)
mean(na.omit(extract(coldestQuart,Eng)[[1]]))

extract(warmestQuart,ACT,fun=mean)
mean(na.omit(extract(warmestQuart,Eng)[[1]]))


# BIO10 = Mean Temperature of Warmest Quarter
# BIO11 = Mean Temperature of Coldest Quarter
# BIO4 = Temperature Seasonality (standard deviation *100)
# BIO7 = Temperature Annual Range (BIO5-BIO6)
# 

