rm(list = ls())
library(raster)
library(circular)
library(doBy)
library(lme4)
library(visreg)
library(reporttools)
library(car)
# #library(CircStats)
# library(gtools)
library(Hmisc)
# 
 #library(plyr)
# library(Hmisc)
library(multcomp)
 library(gplots)

obs.dir<-'/Users/daisy/Google Drive/PhD/Data/Observaitons/Cleaned/Breeding/'

# Breeding quantiles
BP<-read.csv("/Users/daisy/Google Drive/PhD/BreedingTiming/tables/BreedingQuantiles2015-10-17.csv")
# traits
inc<-read.csv('/Users/daisy/Google Drive/PhD/BreedingTiming/tables/SpeciesOfInterest_2015-10-07.csv')
taxon<-inc[c("genus","family","ScientificName")]
BP<-merge(BP,taxon,all.x=TRUE,by.x="Species",by.y="ScientificName")
#peak breeding period data
peak<-read.csv(paste0(obs.dir,"PeakBreedingPointOfLayDayOfYear2015-10-17.csv"))
#make sure only have species interested in
inc2<-subset(inc,remove==1)
RMspecies<-inc2$ScientificName
peak<-peak[peak$Scientific.Name %nin% RMspecies,] #remove unwanted species
BP<-BP[BP$Species %nin% RMspecies,]
#reclassify koeppen zones
koeppen<-raster(paste0('/Users/daisy/Google Drive/PhD/Data/Spatial/BOM_climate_zones/kpngrp_major/koepenReclassified.asc'),
                proj4string=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))
#combine equatorial and tropical: 41 Equatorial, 35 Tropical
m <- c(35, 41, 35)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
koeppen<- reclassify(koeppen,rclmat)
#add koeppen zones to peak data
peak$koeppen<-extract(koeppen,data.frame(cbind(peak$lon, peak$lat)))
#convert data to radians to make circular vector
peak$Radians <-(peak$DOY_PL/366*360)*pi / 180


#=====================
#How does median, start, end, breeding length vary between biomes? 
#=====================

regions<-as.vector(unique(BP$Region))
BSummary<-as.data.frame(paste(regions))
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
#write.csv(RegionSummaryFinal,paste0('/Users/daisy/Google Drive/PhD/BreedingTiming/tables/RegionalSummary',as.Date(Sys.Date()),'.csv'),row.names=FALSE)


#=====================
#ANOVA of 5th 50th and 95th percentils by biome and PBP
#=====================

#adjust dates so that everything less than start is added to 365, 
#then do not need to worry about 1 being equally as close to 365 as 364
BP$MODIFIED_median<-with(BP,ifelse(Quantile50<Quantile5,Quantile50+365,Quantile50))
BP$MODIFIED_end<-with(BP,ifelse(Quantile95<Quantile5,Quantile95+365,Quantile95))


levels(BP$Region) <- abbreviate(levels(BP$Region))

################################
#PBP
################################
  # The syntax of the function call above goes like this:  
  #lmer(outcomeVariable ~ fixedEffect + (1 + fixedEffect | groupingFactorA) + (1 + fixedEffect | groupingFactorB))
#two different lmer with random effects of a) family, b) family and order
modBP1<-lmer(BreedingPeriod ~ Region + (1+Region|family),REML=FALSE,data=BP)#add family
qqPlot(residuals(modBP1),main=("mod1"))#PROFORMS THE BEST BUT NOT GREAT RESIDUALS
modBP2<-lmer(BreedingPeriod ~ Region + (1+Region|family)+(1+Region|Order),REML=FALSE,data=BP)#add family and order, significant but order makes it worse
qqPlot(residuals(modBP2),main=("mod2"))
anova(modBP1,modBP2) #family and order proforms better (loglik)
modANOVA<-Anova(modBP1)#best way to get p values from lmer
#get p, f and Df.res using REML = T
modBP1<-lmer(BreedingPeriod ~ Region + (1+Region|family),REML=T,data=BP)#add family
#formatPval(modANOVA$'Pr(>Chisq)',3,eps=0.001)
Anova(modBP1, test.statistic="F") 
#find differencs between regions
TukeyRegion2<- glht(modBP1, linfct=mcp(Region="Tukey"))
summary(TukeyRegion2)
mn<-summaryBy(BreedingPeriod ~ Region, data=BP,
                FUN=c(mean,SE,sd)) #find mean and sd
lets <- cld(TukeyRegion2)$mcletters$Letters
b <- with(mn, barplot2(BreedingPeriod.mean,
names.arg=Region,
plot.ci=TRUE, ylim=c(0,366),main = "PBP",
ci.l=BreedingPeriod.mean - BreedingPeriod.SE,
ci.u=BreedingPeriod.mean + BreedingPeriod.SE))
text(b, 35, lets)




###################################################

#MEDIAN DATE
#mixed effect model, biome as fixed effect and family is random, 
# p< 0.001
###################################################
BP$Region <- with(BP, reorder(Region,
                              MODIFIED_median, mean))
          
mod1<-lmer(MODIFIED_median ~ Region + (1+Region|family),REML=T,data=BP)#add family
qqPlot(residuals(mod1),main=("mod1"))#PROFORMS THE BEST BUT NOT GREAT RESIDUALS
Anova(mod1, test.statistic="F") 


#mod2<-lmer(MODIFIED_median ~ Region + (1+Region|family)+(1+Region|Order),
 #     REML=FALSE,data=BP)#add family and order, NOT SIGNIFICANT
#qqPlot(residuals(mod2),main=("mod2"))
 #        anova(mod1,mod2)#smaller AIC & BIC the better, make sure P value is significant

# glm1 <- glmer(MODIFIED_median ~ Region + (1+Region|family),data=BP, 
#               family=poisson(link = "identity"))
# qqPlot(residuals(glm1),main=("glm1"))
# nullMod <- glmer(MODIFIED_median ~ (1|family) + (1|Order), 
#                  data=BP,family = poisson)
        # qqPlot(residuals(nullMod),main=("nullMod"))
        # glm2 <- glmer(MODIFIED_median ~ Region + (1+Region|family)+ (1+Region|Order),data=BP,
        #               family = poisson(link = "identity"))#fails to converge
       
modANOVA<-Anova(mod1) #best way to get p values from lmer
modANOVA$'Pr(>Chisq)'
formatPval(modANOVA$'Pr(>Chisq)',3,eps=0.001)

TukeyRegion2<- glht(mod1, linfct=mcp(Region="Tukey"))
lets <- cld(TukeyRegion2)$mcletters$Letters
SE <- function(x)sd(x)/sqrt(length(x))
mn <- summaryBy(MODIFIED_median ~ Region, data=BP,
                FUN=c(mean,SE))
par(mar=c(10,4,4,4), cex.lab=0.7, las=2, cex.axis=0.7)
b <- with(mn, barplot2(MODIFIED_median.mean,
                       names.arg=Region,
                       plot.ci=TRUE, ylim=c(0,366),main="median",
                       ci.l=MODIFIED_median.mean - MODIFIED_median.SE,
                       ci.u=MODIFIED_median.mean + MODIFIED_median.SE))
text(b, 35, lets)
#check on significant levels for desert region and grassland
summary(TukeyRegion2)

mn <- summaryBy(MODIFIED_median ~ Region, data=BP,
                FUN=c(mean,SE))


###################################################

#PBP end
#desert and grassland end PBP significantly earlier than other biomes, 
#desert significantly earlier than other grassland

####################################################
# lm3<-lm(formula = MODIFIED_end ~ Region, data = BP)
# #need to add family
# qqPlot(residuals(lm3),main=("lm3"))
mod3<-lmer(MODIFIED_end  ~ Region + (1+Region|family),REML=FALSE,data=BP)#add family
qqPlot(residuals(mod3),main=("mod1"))#PROFORMS THE BEST BUT NOT GREAT RESIDUALS
mod4<-lmer(MODIFIED_end  ~ Region + (1+Region|family)+(1+Region|Order),
           REML=FALSE,data=BP)#add family and order, NOT SIGNIFICANT
qqPlot(residuals(mod4),main=("mod2"))
anova(mod3,mod4)



# glm1 <- glm(MODIFIED_end ~ Region,data=BP, 
#               family=poisson(link = "identity"))
# summary(glm1)
# qqPlot(residuals(glm1),main=("glm1"))
# 
#failed to converge
# glmer1 <- glmer(MODIFIED_end ~ Region + (1+Region|family),data=BP, 
#               family=poisson(link = "identity"))#much better residuals, notice the scale
# 
# glmer2 <- glmer(MODIFIED_end ~ Region + (1+Region|family) + (1+Region|Order),
#                 data=BP, 
#                 family=poisson(link = "identity"))

Anova(mod4)
TukeyRegion3<- glht(mod4, linfct=mcp(Region="Tukey"))
SE <- function(x)sd(x)/sqrt(length(x))

summary(TukeyRegion3)

lets <- cld(TukeyRegion3)$mcletters$Letters
SE <- function(x)sd(x)/sqrt(length(x))
mn <- summaryBy(MODIFIED_end ~ Region, data=BP,
                FUN=c(mean,SE))
par(mar=c(10,4,4,4), cex.lab=0.7, las=2, cex.axis=0.7)

b <- with(mn, barplot2(MODIFIED_end.mean,
                       names.arg=Region,
                       plot.ci=TRUE, main="end",ylim=c(0,366),
                       ci.l=MODIFIED_end.mean - MODIFIED_end.SE,
                       ci.u=MODIFIED_end.mean + MODIFIED_end.SE))
text(b, 35, lets)


##################################################
#start date

##################################################
# lm4<-lm(formula = Quantile5 ~ Region, data = BP) #original model
# qqPlot(lm(formula = Quantile5 ~ Region, data = BP)) #data not quite normal
# hist(BP$Quantile5,breaks=20)
# fligner.test(Quantile5 ~ Region, data=BP)#regions hace different amoutns of variance
# fligner.test(Quantile5 ~ family, data=BP)
#******this is the best model for this data*******
lm4<-lmer(formula = Quantile5 ~ Region+ (1+Region|family), data = BP,REML=T)
qqPlot(residuals(lm4))
lm5<-lmer(formula = Quantile5 ~ Region+ (1+Region|family)+ (1+Region|Order), data = BP)
qqPlot(residuals(lm5))
anova(lm4,lm5) 



Anova(lm4, test.statistic="F") 





# no significant differnt, use only family
# # #glm with poisson distribution
# glm1<-(glm(Quantile5 ~ Region,data=BP,family = poisson(link = "identity"))) # residuals 14571
# # qqPlot(residuals(glm1))
# # summary(glm1) # 14571

# #glm with poisson distribution, "1" gets (potentially) different intercept for each group in A,
# #fixedEffect" part within the random effect gets you a (potentially) different
# # slope for each group in A. outcomeVariable ~ fixedEffect + (1 + fixedEffect | groupingFactorA)
# glme2 <- glmer(Quantile5 ~ Region + (1+Region|family),data=BP,
#               family = poisson(link = "identity")) 

#find out which model proforms best #LM$
# anova(lm4,glm1, glm2)
summary(lm4)
qqPlot(residuals(lm4))
Anova(lm4)
TukeyRegion4<- glht(lm4, linfct=mcp(Region="Tukey"))
mn <- summaryBy(Quantile5 ~ Region, data=BP,
                FUN=c(mean,SE))
summary(TukeyRegion4)

lets <- cld(TukeyRegion4)$mcletters$Letters
SE <- function(x)sd(x)/sqrt(length(x))
mn <- summaryBy(Quantile5 ~ Region, data=BP,
                FUN=c(mean,SE))
par(mar=c(10,4,4,4), cex.lab=0.7, las=2, cex.axis=0.7)

b <- with(mn, barplot2(Quantile5.mean,
                       names.arg=Region,
                       plot.ci=TRUE, main="start",ylim=c(0,366),
                       ci.l=Quantile5.mean - Quantile5.SE,
                       ci.u=Quantile5.mean + Quantile5.SE))
text(b, 35, lets)



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





