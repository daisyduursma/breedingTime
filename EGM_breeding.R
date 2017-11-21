library(lme4)


#Equally good days of the months

  #=====================================
#funciton for EGM
 
CalcEGM <- function(monthlyObs){
    a<-monthlyObs/sum(monthlyObs)
    b<- log(a)
    b[is.infinite(b)] <- 0 
    EGM<-exp( -( sum(a*b) ) )
    return(EGM)
  }




#==================================
 #peak data
peakSubset <- read.csv("//ad.uws.edu.au/dfshare/HomesHWK$/30022860/Desktop/Daisy/CB/PeakBreedingPointOfLayDayOfYear2015-10-17.csv")
a<-unique(paste(peakSubset$Scientific.Name,peakSubset$koeppen))

peak <- read.csv("//ad.uws.edu.au/dfshare/HomesHWK$/30022860/Desktop/Daisy/CB/PointOfLayDayOfYear2015-10-07.csv")
peak$namepeak<-(paste(peak$Scientific.Name,peak$koeppen))


koeppen<-raster("//ad.uws.edu.au/dfshare/HomesHWK$/30022860/Desktop/Daisy/CB/koepenReclassified.asc",
                proj4string=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))
#combine equatorial and tropical: 41 Equatorial, 35 Tropical
m <- c(35, 41, 35)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
koeppen<- reclassify(koeppen,rclmat)

peak$koeppen<-extract(koeppen,data.frame(cbind(peak$lon, peak$lat)))

peak$namepeak<-(paste(peak$Scientific.Name,peak$koeppen))
peak<-peak[peak$namepeak %in% a,]


#empty list
EGM_ALL<-list()

#for each species
for ( i in 1: length(unique(peak$Scientific.Name))){
  
  sp<-subset(peak,Scientific.Name==unique(peak$Scientific.Name)[[i]]) 
  EGM_sp<-list()
  for(ii in 1:length(unique(sp$koeppen))){
    spbiome<-subset(sp,koeppen==unique(sp$koeppen)[[ii]])
    obsBiome<-as.data.frame(table(spbiome$month))$Freq
    obsBiome<- if (length(obsBiome)==12) { obsBiome
    } else c(obsBiome,rep(0,12-length(obsBiome)))
    EGM_spBiome<-as.data.frame(ifelse(nrow(spbiome)>=40,CalcEGM(obsBiome),NA))
    EGM_spBiome$Species<-unique(peak$Scientific.Name)[[i]]
    EGM_spBiome$Biome<-unique(sp$koeppen)[[ii]]
    colnames(EGM_spBiome)<-c("EGM","Species","Biome")
    EGM_sp[[ii]]<-EGM_spBiome
  }
  
  EGM_ALL[[i]]<-do.call("rbind",EGM_sp)
  
}  

EGMdf<-do.call("rbind",EGM_ALL)
  
EGMdf<-subset(EGMdf, !is.na(EGMdf$EGM))





#Koeppen zones

#   35 Tropical and Equatorial
#   32 Subtropical
#   22 Desert
#   13 Grassland
#   3 Temperate


##############
#Are the differences significant?
##############

#add family to dataset
traits<-read.csv("//ad.uws.edu.au/dfshare/HomesHWK$/30022860/Desktop/Daisy/CB/Australian_Bird_Data_Version_1_0.csv")

traits<-traits[,c("Taxon.scientific.name","Taxon.common.name","Family.scientific.name")]

new<-merge(EGMdf,traits,by.x="Species",by.y="Taxon.scientific.name",all.x=TRUE)

new<-edit(new)

aggregate(EGMdf,by=list(EGMdf$Biome), FUN="mean")# outcome~gender+group, experiment, mean )
aggregate(EGMdf,by=list(EGMdf$Biome), FUN="sd")
levels(EGMdf$Biome) <- abbreviate(levels(EGMdf$Biome))
          

# The syntax of the function call above goes like this:  
#lmer(outcomeVariable ~ fixedEffect + (1 + fixedEffect | groupingFactorA) + (1 + fixedEffect | groupingFactorB))
#two different lmer with random effects of a) family, b) family and order


library(car)
library(multcomp)

new$Biome<-as.factor(new$Biome)
levels(new$Biome) <- abbreviate(levels(new$Biome))

mod1<-lmer(EGM ~ Biome + (1+Biome|Family.scientific.name),REML=F,data=new)
Anova(mod1)#best way to get p values from lmer

# mod2<-lmer(EGM ~ Biome + (1+Biome|Family.scientific.name),REML=T,data=new)#add family
# Anova(mod2, test.statistic="F") 

# qqPlot(residuals(mod1),main=("mod1"))


#get p, f and Df.res using REML = T

TukeyRegion2<- glht(mod1, linfct=mcp(Biome="Tukey"))

summary(TukeyRegion2)

###RESULT - the desert region's EGM is not significantly less or more constrained than any region except tropics and subtropics

write.csv(new,paste0("//ad.uws.edu.au/dfshare/HomesHWK$/30022860/Desktop/Daisy/CB/EGM_breeding",as.Date(Sys.Date()),'.csv'),row.names=FALSE)


################################################
#limit to just species that occur in the desert
################################################
# 
# 
# desSP<-subset(new,Biome==22)$Species
# sharedSP<-new[new$Species %in% desSP,]
# 
# sharedSP$Biome<-as.factor(sharedSP$Biome)
# levels(sharedSP$Biome) <- abbreviate(levels(sharedSP$Biome))
# 
# mod1<-lmer(EGM ~ Biome + (1+Biome|Family.scientific.name),REML=F,data=sharedSP)
# Anova(mod1)#best way to get p values from lmer
# 
# mod2<-lmer(EGM ~ Biome + (1+Biome|Family.scientific.name),REML=T,data=sharedSP)#add family
# Anova(mod2, test.statistic="F") 
# 
# qqPlot(residuals(mod2),main=("mod1"))
# 
# Anova(mod2)#best way to get p values from lmer
# 
# #get p, f and Df.res using REML = T
# 
# TukeyRegion2<- glht(mod2, linfct=mcp(Biome="Tukey"))
# 
# summary(TukeyRegion2)
# 
# 
# 
# 
# 

#


                                                                                  value, mean))


lmClades<-lm(formula = value ~ variable, data = finEGM2)
anova(lmClades)
TukeyRegion2<- glht(lmClades, linfct=mcp(variable="Tukey"))
lets <- cld(TukeyRegion2)$mcletters$Letters
library(doBy)
SE <- function(x)sd(x)/sqrt(length(x))
mn <- summaryBy(value ~ variable, data=finEGM2,
                                 FUN=c(mean,SE))
par(mar=c(10,4,4,4), cex.lab=0.7, las=2, cex.axis=0.7)
b <- with(mn, barplot2(value.mean,
                                               names.arg=variable,
                                               plot.ci=TRUE, ylim=c(0,9),
                                               ci.l=value.mean - value.SE,
                                               ci.u=value.mean + value.SE))
text(b, 8, lets)

  
  write.csv(finEGM,paste0("/Users/daisy/Google Drive/PhD/BreedingTiming/tables/EGM_breeding",as.Date(Sys.Date()),'.csv'),row.names=FALSE)