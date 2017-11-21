
rm(list = ls())
library(raster)
library(data.table)
library(Hmisc)
library(scatterplot3d)
library(lme4)
library(car)
library(multcomp)
library(reporttools)
library(doBy)

obs.dir<-'/Users/daisy/Google Drive/PhD/Data/Observaitons/Cleaned/Breeding/'

# # Breeding quantiles
# BP<-read.csv("/Users/daisy/Google Drive/PhD/BreedingTiming/tables/BreedingQuantiles2015-10-17.csv")
# traits
inc<-fread('/Users/daisy/Google Drive/PhD/BreedingTiming/tables/SpeciesOfInterest_2015-10-07.csv',data.table=FALSE)
#peak breeding period data
peak<-fread(paste0(obs.dir,"PeakBreedingPointOfLayDayOfYear2015-10-17.csv"),data.table=FALSE)
#make sure only have species interested in
inc2<-subset(inc,remove==1)
RMspecies<-inc2$ScientificName
peak<-peak[peak$Scientific.Name %nin% RMspecies,] #remove unwanted species

inc<-inc[c("ScientificName","Jetz.Patch.clade","Common.me","genus","family","order")]

#fix day and month
peak$day<-as.numeric(strftime(as.Date(peak$DOY_PL, origin = "1970-01-01"),format = "%d"))
peak$month<-as.numeric(strftime(as.Date(peak$DOY_PL, origin = "1970-01-01"),format = "%m"))
peak<-merge(peak,inc, by.x= "Scientific.Name",by.y="ScientificName")
#peak$day<-formatC(peak$day,width=2,flag='0')
#peak$month<-formatC(peak$month,width=2,flag='0')

##reclassify koeppen zones
koeppen<-raster(paste0('/Users/daisy/Google Drive/PhD/Data/Spatial/BOM_climate_zones/kpngrp_major/koepenReclassified.asc'),
                proj4string=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))
#combine equatorial and tropical: 41 Equatorial, 35 Tropical
m <- c(35, 41, 35)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
koeppen<- reclassify(koeppen,rclmat)
#add koeppen zones to peak data
peak$koeppen<-extract(koeppen,data.frame(cbind(peak$lon, peak$lat)))
#convert data to radians to make circular vector
#peak$Radians <-(peak$DOY_PL/366*360)*pi / 180

#combine climate and peak obs
clim<-fread(paste0("/Users/daisy/Google Drive/PhD/BreedingTiming/tables/",
                   "SiloClimateLocationDate30DaySummary2015-11-11.csv"),
            data.table = FALSE)
colnames(clim)[2]<-"meanTmax"
MClim<-fread(paste0("/Users/daisy/Google Drive/PhD/BreedingTiming/tables/",
                    "SiloClimateLocationAllMonths2015-11-03.csv"),
             data.table = FALSE)

#add koeppen zones to MCLIM
locs<-as.data.frame(unique(cbind(MClim$lon,MClim$lat)))
locs$koeppen<-extract(koeppen,locs)
MClim2<-merge(MClim,locs,by.x=c("lon","lat"),by.y=c("V1","V2"))
MClim2<-na.omit(MClim2)
#remove white space
clim$day<-as.vector(apply(as.data.frame(clim[,"day"]),2,function(x)gsub('\\s+', '',x)))
clim$month<-as.vector(apply(as.data.frame(clim[,"month"]),2,function(x)gsub('\\s+', '',x)))
clim$ID<-with(clim,paste(ID,day,month,year))
peak$ID<-with(peak,(paste(lat,lon,day,month,year)))
peak<-subset(peak,year>=1957&year <=2013)
climPK<-merge(peak,clim, by.x = "ID",by.y="ID")


#############
#averages by region (control by region)
#############
SE <- function(x)sd(x)/sqrt(length(x))

#precip for species/biome
pr <- summaryBy(sumPrec ~ koeppen + Scientific.Name, data=climPK,
                FUN=c(mean,SE,sd))
#precip by biome
pr<- summaryBy(sumPrec.mean ~ koeppen, data=pr,
               FUN=c(mean,SE,sd))
mnT <- summaryBy(meanTmin ~ koeppen + Scientific.Name, data=climPK,
                 FUN=c(mean,SE,sd))
mnT<- summaryBy(meanTmin.mean ~ koeppen, data=mnT,
                FUN=c(mean,SE,sd))
mxT <- summaryBy(meanTmax ~ koeppen + Scientific.Name, data=climPK,
                 FUN=c(mean))
mxT<- summaryBy(meanTmax.mean ~ koeppen, data=mxT,
                FUN=c(mean,SE,sd))




###############
#Are there differences in temperatures and precipitation, 
#between biomes when accounting for species and family as random effects?
###############
#library(visreg)

climPK$koeppen <- factor(climPK$koeppen)

#TMIN
#mod<-lmer(meanTmin ~ koeppen + (1+koeppen|Scientific.Name),data=climPK)
mod1<-lmer(meanTmin ~ koeppen + (1+koeppen|Scientific.Name)+ (1+koeppen|family),data=climPK, REML=T)
Anova(mod1,test.statistic="F")
#mod2<-lmer(meanTmin ~ koeppen + (1+koeppen|Scientific.Name)+ (1+koeppen|family)+ (1+koeppen|order),data=climPK)
#anova(mod,mod1,mod2,mod3) #compare models

summary(mod1)
modTminANOVA<-Anova(mod1) #best way to get p values from lmer
modANOVA$'Pr(>Chisq)'
formatPval(modANOVA$'Pr(>Chisq)',3,eps=0.005)
TukeyTmin<-glht(mod1, linfct=mcp(koeppen="Tukey"))

#TMax
mod2<-lmer(meanTmax ~ koeppen + (1+koeppen|Scientific.Name)+ (1+koeppen|family),data=climPK)
summary(mod2)
#calculate R2 following = Johnson, Paul C.D. 2014. Extension of Nakagawa & Schielzeth’s R2GLMM to random slopes models. Methods in Ecology and Evolution. DOI: 10.1111/2041-210X.12225.
#rsquared.glmm(list(mod)) #proportion of variance explained by koeppen zone and random effects (temperature) (marginal,conditional)
modTmaxANOVA2<-Anova(mod2) #best way to get p values from lmer
modANOVA2$'Pr(>Chisq)'
formatPval(modANOVA2$'Pr(>Chisq)',3,eps=0.005)
TukeyTmax<-glht(mod2, linfct=mcp(koeppen="Tukey"))


#Precip
mod3<-lmer(sumPrec ~ koeppen + (1+koeppen|Scientific.Name),data=climPK)
summary(mod3)
#calculate R2 following = Johnson, Paul C.D. 2014. Extension of Nakagawa & Schielzeth’s R2GLMM to random slopes models. Methods in Ecology and Evolution. DOI: 10.1111/2041-210X.12225.
#rsquared.glmm(list(mod)) #proportion of variance explained by koeppen zone and random effects (temperature) (marginal,conditional)
modANOVA3<-Anova(mod3) #best way to get p values from lmer
modANOVA3$'Pr(>Chisq)'
formatPval(modANOVA3$'Pr(>Chisq)',3,eps=0.005)
TukeyPrec<-glht(mod3, linfct=mcp(koeppen="Tukey"))







##############
# 5 pannel, width of single-column
#figures should be no less than 1,040 pixels wide
##############

#   35 Tropical
#   32 Subtropical
#   22 Desert
#   13 Grassland
#   3 Temperate

df<-climPK[,c("koeppen","meanTmax","sumPrec","meanTmin")]
colnames(df)[2]<-"maxTemp"

#function to make coloramp transparent     
colorRampAlpha <- function(..., n, alpha) {
  colors <- colorRampPalette(...)(n)
  paste(colors, sprintf("%x", ceiling(255*alpha)), sep="")
}

# Simple function for placing labels on a figure.
plotlabel <- function(txt, where, inset=0.15, font=1, inset.x=inset, inset.y=inset,...){
  u <- par()$usr
  if(grepl("left",where))x <- u[1] + inset.x*(u[2]-u[1])
  if(grepl("right",where))x <- u[2] - inset.x*(u[2]-u[1])
  if(grepl("bottom",where))y <- u[3] + inset.y*(u[4]-u[3])
  if(grepl("top",where))y <- u[4] - inset.y*(u[4]-u[3])
  
  text(x,y,txt,font=font,...)
}




plotNiche<-function(dfFUN,var1,var2,koep,obs,obsVar1,obsVar2,col1,col2,region){
  dfFUN<-subset(dfFUN,koeppen==koep)
#   #####TMIN####
#   x <- densCols(dfFUN[,paste(var1)],dfFUN[,"minTemp"], 
#                 colramp=colorRampPalette(c("black", "white")))
#   dfFUN$dens <- col2rgb(x)[1,] + 1L
#   dfFUN2<-dfFUN[order(- dfFUN$dens,decreasing = FALSE), , drop = FALSE]                                   
#   dfFUN2<-dfFUN2[c(1:(nrow(dfFUN2)*.75)),]#3/4 of observations in denses location
  #get mean of potential and realized
  pTMAX<-mean(dfFUN[,paste(var2)])
  rTMAX<-mean(obs[,paste(obsVar2)])
 
#   ####TMAX
#   #extract 3/4 observations 
#   x <- densCols(dfFUN[,paste(var1)],dfFUN[,paste(var2)], 
#                 colramp=colorRampPalette(c("black", "white")))
#   dfFUN$dens <- col2rgb(x)[1,] + 1L
#   dfFUN2<-dfFUN[order(- dfFUN$dens,decreasing = FALSE), , drop = FALSE]                                   
#   dfFUN2<-dfFUN2[c(1:(nrow(dfFUN2)*.75)),]#3/4 of observations in denses location
#   #get mean 3/4of observations
#   potentialRAIN<-mean(dfFUN2[,paste(var1)])
#   potentialTMAX<-mean( dfFUN2[,paste(var2)])
  #make smaller data sets
  dfFUN[,paste(var1)]<-round(dfFUN[,paste(var1)],digits=1)
  dfFUN[,paste(var2)]<-round(dfFUN[,paste(var2)],digits=1)
  dfFUN<-unique(dfFUN[,c(paste(var1),paste(var2))])
  #plot climate niche
  t.plot<-plot((dfFUN[,paste(var1)]+1),dfFUN[,paste(var2)],
               col='grey80',
               pch=16,
               cex=1,
               ylab="",
               xlab = "",
               ylim=c(-10,50),
               xlim = c(0,2000),
               axes=FALSE
              )
  Axis(side=1, labels=FALSE)
  Axis(side=2, at=c(-10,0,10,20,30,40,50),labels=c(-10,0,10,20,30,40,50))
#   #get TMIN values
#   #extract 3/4 actual observations 
#   x <- densCols(obs[,paste(obsVar1)],obs[,"meanTmin"], 
#                 colramp=colorRampPalette(c("black", "white")))
#   obs$dens <- col2rgb(x)[1,] + 1L
#   obs2<-obs[order(- obs$dens,decreasing = FALSE), , drop = FALSE]                                   
#   obs2<-obs2[c(1:(nrow(obs2)*.75)),]#3/4 of observations in denses location
#   #get mean 3/4of observations
#   actualTMINRAIN<-mean(obs2[,paste(obsVar1)])
#   actualTMIN<-mean( obs2[,"meanTmin"])
  #extract 3/4 actual observations 
  x <- densCols(obs[,paste(obsVar1)],obs[,paste(obsVar2)], 
                                colramp=colorRampPalette(c("black", "white")))
  obs$dens <- col2rgb(x)[1,] + 1L
  obs2<-obs[order(- obs$dens,decreasing = FALSE), , drop = FALSE]                                   
  obs2<-obs2[c(1:(nrow(obs2)*.75)),]#3/4 of observations in denses location

  #add realized to graph
  points((obs[,paste(obsVar1)]+1),obs[,paste(obsVar2)],pch=16,cex=1, col=col1)
  points((obs2[,paste(obsVar1)]+1),obs2[,paste(obsVar2)],pch=16,cex=1, col=col2)
  
  #add lines for means of tmax
  rTMAX<-mean(obs[,paste(obsVar2)])
  abline(h = rTMAX,  col = "black",lwd = 2)
  abline(h = pTMAX, col = "black", lwd = 2,lty = 3)
  #redo x,y axis
  abline(v = 0, col = "black")
  abline(v =-10, col = "black")
  return(t.plot)
}




pdf(file = "/Users/daisy/Google Drive/PhD/BreedingTiming/figures/FEDandTemperature20160127.pdf",
     width = 3.5, height = 9.7)



par(mfrow = c(5,1),
          oma = c(5,4,0,0) + 0.1,
          mar = c(1,2,1,1) + 0.1,
          xaxs="i",
          yaxs = "i", 
          cex=.7)


#Koeppen zones: 35-Tropical, 32-Subtropical, 22-Desert,13-Grassland,3-Temperate
#order of koeppen classes match colours
koepClass<-c(3,13,22,35,32)
kpColour<-c("#1d4e89",
            '#f69256',
            '#401d00',
            '#339900',
            '#7dd0b6')
# ######## Temperate region  
tempSum<-plotNiche(var1="pre",
          var2="maxTemp",
          koep=3, 
          dfFUN=MClim2,
          obs=subset(df,koeppen==3),
          obsVar1="sumPrec",
          obsVar2="maxTemp",
          col1="#33CCFF",
          col2="#1d4e89",
          region="Temperate")
plotlabel("Temperate", "topright")

######## Desert Region
desertSum<-plotNiche(var1="pre",
          var2="maxTemp",
          koep=22, 
          dfFUN=MClim2,
          obs=subset(df,koeppen==22),
          obsVar1="sumPrec",
          obsVar2="maxTemp",
          col1="#D0977B",
          col2="#401d00",
          region="desert")

plotlabel("Desert", "topright")

########## Grasslands
grassSum<-plotNiche(var1="pre",
          var2="maxTemp",
          koep=13, 
          dfFUN=MClim2,
          obs=subset(df,koeppen==13),
          obsVar1="sumPrec",
          obsVar2="maxTemp",
          col1="#FFCC99",
          col2="#f69256",
          region="grassland")

plotlabel("Grassland", "topright")
########## subTropics
subTropSum<-plotNiche(var1="pre",
          var2="maxTemp",
          koep=32, 
          dfFUN=MClim2,
          obs=subset(df,koeppen==32),
          obsVar1="sumPrec",
          obsVar2="maxTemp",
          col1="#CCFFFF",
          col2="#7dd0b6",
          region="Subtropics")
plotlabel("Subtropical", "topright")

########## Tropics
tropicSum<-plotNiche(var1="pre",
          var2="maxTemp",
          koep=35, 
          dfFUN=MClim2,
          obs=subset(df,koeppen==35),
          obsVar1="sumPrec",
          obsVar2="maxTemp",
          col1="#95FF4F",
          col2="#339900",
          region="tropics")
Axis(side=1, at=c(0,500,1000,1500,2000),labels=c(0,500,1000,1500,2000))

#add x and y labels
mtext(text = "Precipitation (mm)", 1, 2, outer=TRUE)
mtext(text = expression(paste("Maximum temperature (",degree,"C)")),
      2, 2, outer=TRUE,las=0)


plotlabel("Tropical", "topright")
dev.off()




##########################

# only Desert


##################################
plotNiche<-function(dfFUN,var1,var2,koep,obs,obsVar1,obsVar2,col1,col2,region){
  dfFUN<-subset(dfFUN,koeppen==koep)

  #get mean of potential and realized
  pTMAX<-mean(dfFUN[,paste(var2)])
  rTMAX<-mean(obs[,paste(obsVar2)])
  #make smaller data sets
  dfFUN[,paste(var1)]<-round(dfFUN[,paste(var1)],digits=1)
  dfFUN[,paste(var2)]<-round(dfFUN[,paste(var2)],digits=1)
  dfFUN<-unique(dfFUN[,c(paste(var1),paste(var2))])
  #plot climate niche
  t.plot<-plot((dfFUN[,paste(var1)]+1),dfFUN[,paste(var2)],
               col='grey80',
               pch=16,
               cex=1,
               ylab="",
               xlab = "",
               ylim=c(5,50),
               xlim = c(0,1500),
               axes=FALSE
  )
  Axis(side=2, at=c(10,20,30,40,50),labels=c(10,20,30,40,50))

  #extract 3/4 actual observations 
  x <- densCols(obs[,paste(obsVar1)],obs[,paste(obsVar2)], 
                colramp=colorRampPalette(c("black", "white")))
  obs$dens <- col2rgb(x)[1,] + 1L
  obs2<-obs[order(- obs$dens,decreasing = FALSE), , drop = FALSE]                                   
  obs2<-obs2[c(1:(nrow(obs2)*.75)),]#3/4 of observations in denses location
  
  #add realized to graph
  points((obs[,paste(obsVar1)]+1),obs[,paste(obsVar2)],pch=16,cex=1, col=col1)
  points((obs2[,paste(obsVar1)]+1),obs2[,paste(obsVar2)],pch=16,cex=1, col=col2)
  
  #add lines for means of tmax
  rTMAX<-mean(obs[,paste(obsVar2)])
  abline(h = rTMAX,  col = "black",lwd = 2)
  abline(h = pTMAX, col = "black", lwd = 2,lty = 3)
  #redo x,y axis
  abline(v = 0, col = "black")
  abline(v =-10, col = "black")
  return(t.plot)
}




pdf(file = "/Users/daisy/Google Drive/PhD/BreedingTiming/figures/DesertFEDandTemperature20160127.pdf",
    width = 3, height = 3)



par(mfrow = c(1,1),
    oma = c(5,4,0,0) + 0.1,
    mar = c(1,2,1,1) + 0.1,
    xaxs="i",
    yaxs = "i", 
    cex=.7)


#Koeppen zones: 35-Tropical, 32-Subtropical, 22-Desert,13-Grassland,3-Temperate
#order of koeppen classes match colours
koepClass<-22

######## Desert Region
desertSum<-plotNiche(var1="pre",
                     var2="maxTemp",
                     koep=22, 
                     dfFUN=MClim2,
                     obs=subset(df,koeppen==22),
                     obsVar1="sumPrec",
                     obsVar2="maxTemp",
                     col1="#D0977B",
                     col2="#401d00",
                     region="desert")

Axis(side=1, at=c(0,500,1000,1500),labels=c(0,500,1000,1500))

#add x and y labels
mtext(text = "Precipitation (mm)", 1, 2, outer=TRUE)
mtext(text = expression(paste("Maximum temperature (",degree,"C)")),
      2, 2, outer=TRUE,las=0)

dev.off()




