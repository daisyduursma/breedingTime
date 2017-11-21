rm(list = ls())
library(raster)
library(circular)
library(doBy)
library(lme4)
#library(visreg)
#library(reporttools)
library(car)
# #library(CircStats)
# library(gtools)
library(Hmisc)
# 
 #library(plyr)
# library(Hmisc)
library(multcomp)
#library(gplots)
library(heavy)
obs.dir<-'/Users/daisy/Google Drive/PhD/Data/Observaitons/Cleaned/Breeding/'
# Breeding quantiles
BP<-read.csv(paste0(obs.dir,'PointOfLayDayOfYear2016-09-20.csv'))
# traits
inc<-read.csv('/Users/daisy/Google Drive/PhD/BreedingTiming/tables/SpeciesOfInterest_2016-10-07.csv')
taxon<-inc[c("genus","Clem.Family","ScientificName","Clements.Scientific.name","sort.v2016","family","English.name","Order")]
BP<-merge(BP,taxon,all.x=TRUE,by.x="Scientific.Name",by.y="ScientificName")



########################################################
#loop through each region get breeding observations so preference is given to egg, young, and multi
#######################################################

regions<-unique(BP$koeppen)
BPregion<-list()
for(rg in 1:length(regions)){#loop through differnt koeppen zones
  # #keep only temperate region
  dat<-subset(BP,koeppen==regions[rg])
  #get species with more than 100 Egg or young or multi observations
  datAcc<-subset(dat, type == "solitary-egg" | type == "multi"| type == "young")
  datAcc_sp<-as.data.frame(table(datAcc$Scientific.Name))
  datAcc_sp<-as.vector(droplevels(subset(datAcc_sp,Freq>=100)$Var1))
  datAcc_sp<-data.frame(Species =datAcc_sp, group ="HighAcc")
  datAcc<-datAcc[datAcc$Scientific.Name %in% datAcc_sp$Species,]
  #get species with more than 50  observations and include all data
  dat2<-dat[dat$Scientific.Name %nin% datAcc_sp$Species,]
  dat2_sp<-as.data.frame(table(dat2$Scientific.Name))
  dat2_sp<-as.vector(droplevels(subset(dat2_sp,Freq>=50)$Var1))
  dat2<-dat2[dat2$Scientific.Name %in% dat2_sp,]
  #put data back together
  BPregion[[rg]]<-rbind(datAcc,dat2)
}
obs<-do.call("rbind",BPregion)  

########################################################
#calalculate breeding summaries 
#######################################################
species<-unique(obs$Scientific.Name)
datSummary<-list()
for (i in 1:length(species)){#species data
  spdat<-subset(obs,Scientific.Name==as.character(species[i]))
  #get koeppens for species
  koeps<-unique(spdat$koeppen)
  kpSummary<-list()
  for(kp in 1:length(koeps)){#koepen level species data
    kpdat<-subset(spdat,koeppen ==koeps[kp] )#get species data
    #calculate summary statistics
    NoObs<-as.numeric(nrow(kpdat))#get observation count
    obsRadians<-kpdat$Radians
    #date of quantiles
    quant<-quantile.circular(obsRadians,c(.05,.5,.95),type=8)
    Per5<-round((quant[[1]]*180/pi)/360*365)
    Per50<-round((quant[[2]]*180/pi)/360*365)
    Per95<-round((quant[[3]]*180/pi)/360*365)
    #Breeding Period Length
    StartDate<-as.numeric(Per5)
    EndDate<-as.numeric(Per95)
    BPL<-ifelse (EndDate < StartDate, 365-StartDate+EndDate,EndDate - StartDate )
    #includes unknown
    unKn<-nrow(subset(kpdat,type=="unknown"))
    # #put the data together and name it
    Perdat<-cbind(as.vector(species[i]),koeps[kp],NoObs,
                  Per5,Per50,Per95,BPL,unKn>0)
    colnames(Perdat)<-c('Species','Region','Observations','Quantile5',
                        'Quantile50','Quantile95','BreedingPeriod','dataIncludesUnknown')
    kpSummary[[kp]]<-Perdat
  }

  datSummary[[i]]<-as.data.frame(do.call('rbind', kpSummary))

}

regSummary<-as.data.frame(do.call("rbind",datSummary))

#add in Clements taxonomy and write out

regSummary<-merge(regSummary,taxon,all.x=TRUE,by.x="Species",by.y="ScientificName")
#write.csv(regSummary,"/Users/daisy/Google Drive/PhD/BreedingTiming/manuscript/TablesFigures/SupplementaryData12016_10_11.csv",row.names=FALSE)


# #=====================
# #How does median, start, end, breeding length vary between biomes? 
# #=====================


regSummary$BreedingPeriod2 <-as.numeric(levels(regSummary$BreedingPeriod))[regSummary$BreedingPeriod]
regSummary$Quantile5 <-as.numeric(levels(regSummary$Quantile5))[regSummary$Quantile5]
regSummary$Quantile50 <-as.numeric(levels(regSummary$Quantile50))[regSummary$Quantile50]
#find median for peak
medRadian <-(regSummary$Quantile50/365*360)*pi / 180
medquant<-quantile.circular(medRadian,c(.5),type=8)
medDate<-round((medquant*180/pi)/360*365)
regSummary$Centered50<-ifelse(regSummary$Quantile50<99,regSummary$Quantile50+265,regSummary$Quantile50-100)
#find median for start
stmedRadian <-(regSummary$Quantile5/365*360)*pi / 180
stmedquant<-quantile.circular(stmedRadian,c(.5),type=8)
stmedDate<-round((stmedquant*180/pi)/360*365)
regSummary$Centered5<-ifelse(regSummary$Quantile5<35,regSummary$Quantile5+335,regSummary$Quantile5-36)




regSummary$Quantile95 <-as.numeric(levels(regSummary$Quantile95))[regSummary$Quantile95]



regions<-as.vector(unique(regSummary$Region))
BSummary<-as.data.frame(paste(regions))
colnames(BSummary)<-"Region"
RegionSummary<-list()
for (r in 1:length(regions)){
  D5th<-subset(regSummary,Region==regions[r],select=Quantile5)
  D95th<-subset(regSummary,Region==regions[r],select=Quantile95)
  D50th<-subset(regSummary,Region==regions[r],select=Quantile50)
  vals<-subset(regSummary,Region==regions[r],select=BreedingPeriod)
  birds<-length(vals[,1])#number of breeding birds
  meanPeriod<-mean(as.numeric(vals[,1]))#average length of breeding period
  sdPeriod<-sd(as.numeric(vals[,1])) #sd of breeding perid
  D5thr <-(as.numeric(na.omit(D5th)[,1])/365*360)*pi / 180#median date to radians
  D5thrC<-circular(D5thr,modulo ="2pi", units="radians", rotation="counter") #to circular
  AvgStartR<-mean(circular(D5thr,modulo ="2pi", units="radians", rotation="counter"))
  sdStartR<-sd(circular(D5thr,modulo ="2pi", units="radians", rotation="counter"))
  AvgStart<-round((AvgStartR[[1]]*180/pi)/360*365)#average median date as day
  sdStart<-round((sdStartR[[1]]*180/pi)/360*365)
  D50thr <-(as.numeric(na.omit(D50th)[,1])/365*360)*pi / 180#median date to radians
  D50thrC<-circular(D50thr,modulo ="2pi", units="radians", rotation="counter") #to circular
  AvgMedR<-mean(circular(D50thr,modulo ="2pi", units="radians", rotation="counter"))
  AvgMed<-round((AvgMedR[[1]]*180/pi)/360*365)#average median date as day
  sdMedR<-sd(circular(D50thr,modulo ="2pi", units="radians", rotation="counter"))
  sdMed<-round((sdMedR[[1]]*180/pi)/360*365)
  D95thr <-(as.numeric(na.omit(D95th)[,1])/365*360)*pi / 180#end date to radians
  D95thrC<-circular(D95thr,modulo ="2pi", units="radians", rotation="counter") #to circular
  AvgendR<-mean(circular(D95thr,modulo ="2pi", units="radians", rotation="counter"))
  AvgEnd<-round((AvgendR[[1]]*180/pi)/360*365)#end median date as day
  sdendR<-sd(circular(D95thr,modulo ="2pi", units="radians", rotation="counter"))
  sdEnd<-round((sdendR[[1]]*180/pi)/360*365)#end median date as day


  region<-regions[r]
  RegionSummary[[r]]<-cbind(region,birds,meanPeriod,sdPeriod,AvgStart,sdStart,AvgEnd,sdEnd,AvgMed,sdMed)
}
RegionSummaryFinal<-as.data.frame(do.call("rbind",RegionSummary))#convert to data.frame
# #write.csv(RegionSummaryFinal,paste0('/Users/daisy/Google Drive/PhD/BreedingTiming/tables/RegionalSummary',as.Date(Sys.Date()),'.csv'),row.names=FALSE)



#=====================
#ANOVA of 5th 50th and 95th percentils by biome and PBP
#=====================



#adjust dates so that everything less than start is added to 365, 
#then do not need to worry about 1 being equally as close to 365 as 364
regSummary$MODIFIED_median<-with(regSummary,ifelse(Quantile50<Quantile5,Quantile50+365,Quantile50))
regSummary$MODIFIED_end<-with(regSummary,ifelse(Quantile95<Quantile5,Quantile95+365,Quantile95))



################################
#PBP
################################

regSummary$Region <- with(regSummary, reorder(Region,
                                              BreedingPeriod, mean))

  # The syntax of the function call above goes like this:  
  #lmer(outcomeVariable ~ fixedEffect + (1 + fixedEffect | groupingFactorA) + (1 + fixedEffect | groupingFactorB))

modBP1<-lmer(BreedingPeriod2 ~ Region + (1|Order),data=regSummary)#add family
qqPlot(residuals(modBP1),main=("mod1"))#PROFORMS THE BEST BUT NOT GREAT RESIDUALS
#modANOVA<-Anova(modBP1)#best way to get p values from lmer
#get p, f and Df.res using REML = T
modBP1<-lmer(BreedingPeriod2 ~ Region + (1|Order),REML=T,data=regSummary)#add family
Anova(modBP1, test.statistic="F")
#find differencs between regions
TukeyRegion2<- glht(modBP1, linfct=mcp(Region="Tukey"))
summary(TukeyRegion2)
summaryBy(BreedingPeriod2 ~ Region, data=regSummary,
                FUN=c(mean,sd)) #find mean and sd
 

##################################################
#start date

##################################################

# lmstart<-lmer(formula = Quantile5 ~ Region+ (1+Region|family), data = regSummary,REML=T)
# qqPlot(residuals(lmstart))
# Anova(lmstart, test.statistic="F") 
# TukeyRegion4<- glht(lmstart, linfct=mcp(Region="Tukey"))
# summaryBy(Quantile5 ~ Region, data=regSummary,
#                 FUN=c(mean,sd,length))
# summary(TukeyRegion4)

#----------------------------------
#looking at Orders Start
#----------------------------------
#make new group
regSummary$Region <- factor(as.character(regSummary$Region), levels = c("22", "13","3", "32", "35"))
regSummary$newgroup <- factor(paste0(as.character(regSummary$Region),as.character(regSummary$family)))
hlm<-heavyLme(Centered5 ~ Region, random = ~1, groups = ~newgroup, data=regSummary)
summary(hlm)

OrderdatStart<-list()#empty list
#loop through each order
for(o in 1:length(unique(regSummary$Order))){
  df<-subset(regSummary,Order==unique(regSummary$Order)[o])#get data for order
  if(nrow(df) <2) {#get details if there are not enough observations
    outdat<-t(data.frame(c(as.character(unique(regSummary$Order)[o]),
                           nrow(df),
                           length(unique(df$Species))
    )))
    row.names(outdat) <- NULL 
    colnames(outdat)<-c("Order","numberObservations","numberSpecies")
    OrderdatStart[[o]]<-outdat
  } else {#if there are enough observations find out if significant diffs between biomes
    
    hlm2<-heavyLm(Quantile5 ~ Region, data=df)
    cfDF<-summary(hlm2)$coefficients
    
    #make dates right
   d1<-as.vector(cfDF[,1])
   for(cf in 2:length(d1)){
        d1[cf]<-d1[cf]+d1[1]
   }

    outdat<-t(data.frame(c(as.character(unique(regSummary$Order)[o]),
                           nrow(df),
                           length(unique(df$Species)),
                           paste(round(d1),"(\u00b1",round(cfDF[,2]),")",as.vector(round(cfDF[,4],3)),sep=""),
                           row.names = NULL)))
    row.names(outdat) <- NULL 
    
    colnames(outdat)<-c("Order","numberObservations","numberSpecies",as.character(levels(hlm2$model[,2])))
    
    OrderdatStart[[o]]<-outdat
  }
}

OrderdatStart2<-do.call("smartbind",OrderdatStart)
OrderdatStart2



#----------------------------------
#looking at Orders Peak
#----------------------------------
#make new group
regSummary$Region <- factor(as.character(regSummary$Region), levels = c("22", "13","3", "32", "35"))
regSummary$newgroup <- factor(paste0(as.character(regSummary$Region),as.character(regSummary$family)))
hlm<-heavyLme(MODIFIED_median ~ Region, random = ~1, groups = ~newgroup, data=regSummary)
summary(hlm)

OrderdatPeak<-list()#empty list
#loop through each order
for(o in 1:length(unique(regSummary$Order))){
  df<-subset(regSummary,Order==unique(regSummary$Order)[o])#get data for order
  if(nrow(df) <2) {#get details if there are not enough observations
    outdat<-t(data.frame(c(as.character(unique(regSummary$Order)[o]),
                           nrow(df),
                           length(unique(df$Species))
    )))
    row.names(outdat) <- NULL 
    colnames(outdat)<-c("Order","numberObservations","numberSpecies")
    OrderdatStart[[o]]<-outdat
  } else {#if there are enough observations find out if significant diffs between biomes
    
    hlm2<-heavyLm(MODIFIED_median ~ Region, data=df)
    cfDF<-summary(hlm2)$coefficients
    
    #make dates right
    d1<-as.vector(cfDF[,1])
    for(cf in 2:length(d1)){
      d1[cf]<-d1[cf]+d1[1]
    }
    
    outdat<-t(data.frame(c(as.character(unique(regSummary$Order)[o]),
                           nrow(df),
                           length(unique(df$Species)),
                           paste(round(d1),"(\u00b1",round(cfDF[,2]),")",as.vector(round(cfDF[,4],3)),sep=""),
                           row.names = NULL)))
    row.names(outdat) <- NULL 
    
    colnames(outdat)<-c("Order","numberObservations","numberSpecies",as.character(levels(hlm2$model[,2])))
    
    OrderdatPeak[[o]]<-outdat
  }
}

OrderdatPeak2<-do.call("smartbind",OrderdatPeak)
OrderdatPeak2



#----------------------------------
#looking at families
#----------------------------------


familydat<-list()
for(f in 1:length(unique(regSummary$Clem.Family))){
  D1<-subset(regSummary,Clem.Family==unique(regSummary$Clem.Family)[f])
  if(nrow(D1) <2) {
    next 
  }
  
  if(nrow(D1) <2) {
    next 
  }
  hlm2<-heavyLm(Centered50 ~ Region, data=D1)
  print(paste0(unique(regSummary$Clem.Family)[f],", ",regSummary$Order[1],", number of species = ",nrow(df)))
  print(summary(hlm2))
  
  familydat[[f]]<-cbind(unique(regSummary$Clem.Family)[f],round(summary(hlm2)$coefficients[,4],-3)
  
  # hlm<-heavyLme(Centered50 ~ Region, random = ~1, groups = ~newgroup, data=df)
  # hlm<-heavyLme(Centered50 ~ Region,  random = ~Region, groups = ~family, data=df)
  # summary(hlm)
}

#---------------------------------








  
hlm2<-heavyLm(Centered50 ~ Region, data=regSummary2)
summary(hlm2)

#working
hlm3<-heavyLme(MODIFIED_median ~ Region, random = ~1, groups=~Order, data=regSummary)  


# draws from t distribution
# rand.tdata <- rt(1e4,df=2) #10,000 samples from student t dist with 2 degrees freedom
# 
# # quantile-quantile plot based on normal distribution
# qqnorm(rand.tdata)
# qqline(rand.tdata)
# quantile-quantile plot based on t distribution
# rand.tdata <- sort(rand.tdata)
# plot(qt((1:1e4)/(1e4+1),df=2),rand.tdata)
# abline(0,1)

#residuals not on a t distribution
qqnorm(regSummary$MODIFIED_median)
qqnorm(hlm3$Resid[,2])
qqline(hlm3$Resid[,2])

# res <- sort(hlm3$Resid[,2])
# n <- length(res)
# 
# # weird abline
# plot(qt((1:n)/(n+1),df=2), res)
# abline(0,1)
# 
# 
# # ok
# qqplot(res, rt(n, df=2))
# qqline(res, rt(n, df=2))

# best
qqPlot(hlm3$Resid[,2], distribution="t", df=2)





 









Anova(mod1, test.statistic="F") 


TukeyRegion2<- glht(hlm2, linfct=mcp(Region="Tukey"))

mn <- summaryBy(Centered50 ~ Region, data=regSummary2,
                FUN=c(mean,sd))

mn <- summaryBy(Centered50 ~ Region + family, data=regSummary2,
                FUN=c(mean,sd))

mn2 <- summaryBy(Centered50.mean ~ Region, data=mn,
                FUN=c(mean,sd))

#check on significant levels for desert region and grassland
summary(TukeyRegion2)





# =====================================

# Figure 2. circular plot of Australia with density curves for biomes
# ===================================== 
{
  #Size of figures A single column width measures 88 mm 
  #and a double column width measures 180 mm. 
  #In practice this means that the absolute width of single-column
  #figures should be no less than 1,040 pixels wide
  #and double-column figures should be no less than 2,080 pixels wide 
  #(excluding peripheral white space).
  
# jpeg(file = "/Users/daisy/Google Drive/PhD/BreedingTiming/figures/breedingObservationDensity20151117.jpeg",
#       width = 8, height = 12,units = "cm",res=300)

pdf(file = "/Users/daisy/Google Drive/PhD/BreedingTiming/figures/breedingObservationDensity24151124.pdf",
       width = 3.15, height = 4.73)

  #same par setting as other figures
  par(mfrow=c(2,1), oma=c(0,0,0,0), mar=c(1,1,1,1), las=1, cex=1, xaxs="i",
      yaxs = "i")
  #peak breeding period radions
  PLc<- circular(peak$Radians, units = "radians", zero=pi/2, 
                 rotation = "clock")
  #get midpoints of 12 months on circle
  bins<-12
  arc <- (2 * pi)/bins
  sector <- seq(0, 2 * pi - (2 * pi)/bins, length = bins)
  midpoints <- circular(seq(arc/2, 2 * pi - pi/bins, length = bins), units = "radians", zero=pi/2, 
                        rotation = "clock")
  #rose diagram of Australian breedin species
  rose.diag(PLc, bins=12, col="NA", cex=1.1, prop=1.7,
            axes = FALSE,ticks = FALSE, shrink=1.5,border='NA')
  #add month labels
  axis.circular(at=midpoints, labels=c(1:12), cex=.8,tcl=0.075)
  #Koeppen zones: 35-Tropical, 32-Subtropical, 22-Desert,13-Grassland,3-Temperate
  #order of koeppen classes match colours
  koepClass<-c(3,13,22,35,32)
  kpColour<-c("#1d4e89",
               '#f69256',
               '#401d00',
               '#339900',
               '#7dd0b6')
  #plot kernal density lines for each koeppen region
  lines<-c(1,1,1,1,1)
  for (k in 1:length(koepClass)){
    kpdat<-circular(subset(peak,koeppen==koepClass[k],select = Radians), units = "radians", zero=pi/2, 
                    rotation = "clock")
    lines(density.circular(kpdat, bw=12), lwd=2.5, lty=lines[k],
          col = kpColour[k],shrink=.6)
    #arrows.circular(median(kpdat), lwd=2, lty=lines[k],col = kpColour[k])#arrow pointing to median date
  }
  
  #=======AUS Number of species breeding in any month
  #get unique species per month in Australia and make rose diagram
  monthdat<-peak[!duplicated(peak[c("Scientific.Name","month")]),]
  Mthc<- circular(monthdat$Radians, units = "radians", zero=pi/2, 
                  rotation = "clock")
  rose.diag(Mthc, bins=12, col="NA", cex=1.1, prop=1.7,
            axes = FALSE,ticks = FALSE, shrink=1.5,border='NA')
  axis.circular(at=midpoints, labels=c(1:12), cex=.8,tcl=0.075)
  #for biomes
  for (k in 1:length(koepClass)){
    kpdat<-subset(monthdat,koeppen==koepClass[k])
    kpdat<- kpdat[!duplicated(kpdat[c("Scientific.Name","month")]),]
    kpdat<-circular(kpdat$Radians, units = "radians", zero=pi/2, 
                    rotation = "clock")
    lines(density.circular(kpdat, bw=12), lwd=2.5,
          lty=lines[k],col = kpColour[k],shrink=.6)
    #arrows.circular(median(kpdat), lwd=2, lty=lines[k],col = kpColour[k])
  }
#   text("a",x=-1.5,y=4.5,font=2,cex=1)
#   text("b",x=-1.5,y=1.15, font =2,cex=1)
#   
  
dev.off()
  
  

#biomes with data in

pdf(file = "/Users/daisy/Google Drive/PhD/BreedingTiming/figures/biomes24151117.pdf",
     width = 5.12, height = 4.73)


koeppen2<-reclassify(koeppen, matrix(c(31,33,45),ncol=3))
image (koeppen2,col=kpColour,axes=F,xlab="",ylab="",xlim=c(101,160), ylim=c(-45,-8))  
#temperate
#    kpColour<-c('#ba8c5d',
#               '#339900',
#               '#401d00',
#               '#f4d100',
#               '#186fb8')


kpColour<-c("#1d4e89",
            '#f69256',
            '#401d00',
            '#339900',
            '#7dd0b6')



legend(110,-38,legend="",
       pch=15 ,  bty="n", text.font=1, 
       col = kpColour[3],
       cex=2)
#grassland
legend(131,-38,legend="",
       pch=15 ,  bty="n", text.font=1, 
       col = kpColour[2],
       cex=2)
#desert
legend(151,-38,legend="",
       pch=15 ,  bty="n", text.font=1, 
       col = kpColour[1],
       cex=2)
#tropical
legend(151,-8,legend="",
       pch=15 ,  bty="n", text.font=1, 
       col = kpColour[4],
       cex=2)
#subtropical
legend(151,-23,legend="",
       pch=15 ,  bty="n", text.font=1, 
       col = kpColour[5],
       cex=2)
dev.off()



# ===================================== 

# Table of Biome Summary

# ===================================== 

head(BP)
regions<-unique(BP$Region)
BSummary<-as.data.frame(paste(c("Australia", kpName)))
colnames(BSummary)<-"Region"

RegionSummary<-list()
for (r in 1:length(regions)){
  D5th<-subset(BP,Region==regions[r],select=Quantile5)
  D95th<-subset(BP,Region==regions[r],select=Quantile95)
  D50th<-subset(BP,Region==regions[r],select=Quantile50)
  vals<-subset(BP,Region==regions[r],select=BreedingPeriod)
  birds<-length(vals[,1])#number of breedin birds
  meanPeriod<-mean(as.numeric(vals[,1]))#average length of breeding period
  sdPeriod<-sd(as.numeric(vals[,1])) #sd of breeding perid
  D50thr <-(as.numeric(na.omit(D50th)[,1])/365*360)*pi / 180#median date to radians
  D50thrC<-circular(D50thr,modulo ="2pi", units="radians", rotation="counter") #to circular
  #plot(D50thrC)
  AvgMedR<-mean(circular(D50thr,modulo ="2pi", units="radians", rotation="counter"))
  AvgMed<-round((AvgMedR[[1]]*180/pi)/360*365)#average median date as day
  Rbar<-rho.circular(D50thrC)#average clustering 0 is uncluseted 1 is all at same location
  #message(Rbar)
  V<-1-Rbar#sample circular variance
  #change this to average skewness across all species
  #skewness, need at least 3 data points to calc
  #negative numbers skewed counter clockwise, poitive clockwise, 0 = not skewed
  t2t <- trigonometric.moment(D50thrC, p=2, center=TRUE)
  bbar2 <- t2t$sin
  skewness <- bbar2/(V**(3/2)) #skewness 
  #abar2 <- t2t$cos
  #hatk <- (abar2-Rbar**4)/(V**2) ; #kurtosis
  region<-regions[r]
  RegionSummary[[r]]<-cbind(region,birds,meanPeriod,sdPeriod,AvgMed,skewness,Rbar)
}

RegionSummaryFinal<-as.data.frame(do.call("rbind",RegionSummary))#convert to data.frame
write.csv(RegionSummaryFinal,paste0('/Users/daisy/Google Drive/PhD/BreedingTiming/tables/RegionalSummary',as.Date(Sys.Date()),'.csv'),row.names=FALSE)





