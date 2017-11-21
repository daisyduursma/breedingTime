rm(list = ls())
library(raster)
library(circular)
#library(CircStats)
library(gtools)



dat.dir<-'/Users/daisy/Google Drive/PhD/Data/Observaitons/Cleaned/Breeding/'
#read in observations
finPL<-read.csv(paste0(dat.dir,'PointOfLayDayOfYear2015-07-09.csv'))
#remove species not included in study
#get traits
inc<-read.csv('/Users/daisy/Google Drive/PhD/BreedingTiming/tables/SpeciesOfInterest_2015-07-06.csv')
inc2<-subset(inc,remove==1)
RMspecies<-inc2$ScientificName

finPL<-finPL[finPL$Scientific.Name %nin% RMspecies,]



#make figure of density of observations in years for 3 types of data
NRSYear <-subset(finPL, sourceName == "NestRecordScheme")$year
musYear<-subset(finPL, sourceName == "SouthAustraliaMuseum" |
                  sourceName == "QueenVictoria"| 
                  sourceName ==  "AustraliaMuseum"|
                  sourceName == "MuseumVictoria" |
                  sourceName == "AustralianNationalWildlifeCollection" |
                  sourceName == "WesternAustraliaMuseum"|
                  sourceName == "TasmanianMuseumArtGallery" |
                  sourceName == "NorthernTerritoryMuseumAndArtGallery" |
                  sourceName == "QueenslandMuseum") $year 
altYear<-finPL[grep("^ATLAS", finPL$sourceName), ]$year
ABBBSyear<-subset(finPL, sourceName == "ABBBS")$year
ebirdyear<-subset(finPL, sourceName == "eBIRD")$year
plot(density(altYear), main ="PDF year of collection",xlab="year",xlim=c(1770,2015),col="purple3", lwd=2)
lines(density(NRSYear), col="green", lwd=2)
lines(density(musYear), col="deepskyblue", lwd=2)
lines(density(ABBBSyear), col="darkred", lwd=2)
lines(density(ebirdyear), col="darkgoldenrod2", lwd=2)
legend("topleft",c("Atlas","NRS","Museum","ABBBS","eBIRD"),text.col=c("purple3","green","deepskyblue","darkred","darkgoldenrod2"))

#ATLAS
length(altYear)/nrow(finPL)*100
#Mus
length(musYear)/nrow(finPL)*100
#ebird
length(ebirdyear)/nrow(finPL)*100
#ABBBS
length(ABBBSyear)/nrow(finPL)*100
#NRS
length(NRSYear)/nrow(finPL)*100

###########################

#Peak Breeding Obs

##########################

peak<-read.csv(paste0(dat.dir,"PeakBreedigPointOfLayDayOfYear2015-07-12.csv"))
peak<-peak[peak$Scientific.Name %nin% RMspecies,]
nrow(peak)



#==========================

#temperature and observations by biome

#==========================

library(gtools)

#climate variabls
clim<-c('prec','tmax','tmin' )
climateDir<-"/Users/daisy/Google Drive/PhD/Data/Spatial/Climate/EMASTClimate_mmn"

#koeppen zones
koeppen<-raster(paste0('/Users/daisy/Google Drive/PhD/Data/Spatial/BOM_climate_zones/kpngrp_major/koepenReclassified.asc'),
                proj4string=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))
r1 <- raster("/Users/daisy/Google Drive/PhD/Data/Spatial/Climate/EMASTClimate_mmn/tmin/eMAST_ANUClimate_mmn_tmin_v1m0_01.nc")
#make koeppen zones (more course than climate) so that they have same
#resolution and extent as climate
kp2<-resample(koeppen,r1,method="ngb")
m <- c(35, 41, 35)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
koeppen<- reclassify(kp2, rclmat)


#Peak breeding observations
dat<-read.csv(paste0(dat.dir,"PeakBreedigPointOfLayDayOfYear2015-07-12.csv"))
Locs <- SpatialPoints(cbind(dat$lon, dat$lat),
                      proj4string = CRS('+proj=longlat +datum=WGS84'))

#climatic conditions for biomes
endClimDat<-list()

for (c in 1:length(clim)) {
  #get climate files
  files <- list.files(paste0(climateDir,"/",clim[c]))
  # make a loop to sort through the files
  climdat<-list()
  for (i in 1:length(files)) {
    # make a raster for the desired timeperiod
    r1 <- raster(paste0(climateDir,"/",clim[c],"/",files[i]))
    mr1<-cellStats(r1,mean)
    sdr1<-cellStats(r1,sd)
    month<-strsplit(strsplit(files[i],"_")[[1]][6],".nc")[[1]][1]
    # extract the data and write it directly to the csv file
    MeanBiomes <- zonal(r1,koeppen, fun = mean,na.rm=TRUE)
    sdBiomes <- zonal(r1,koeppen, fun = sd,na.rm=TRUE)
    # reasign the column name to the one we generated
    enddat<-merge(MeanBiomes,sdBiomes, by="zone")
    colnames(enddat)<-c("Biome",paste0(clim[c],"Mean"),paste0(clim[c],"Sd"))
    enddat<-rbind(enddat,c("Australia",mr1,sdr1))
    enddat$Month<-month
    climdat[[i]]<-enddat
  }
  endClimDat[[c]]<-data.frame(do.call("smartbind",climdat))
}
#make table of summary  
final<-do.call("cbind",endClimDat)[c("Biome","Month","precMean","precSd","tmaxMean","tmaxSd",
                                     "tminMean","tminSd")]

#make 2X3 graph of breeding density in biome and precip and temp
par(mfrow=c(2,3), mar = c(3,4,2,3))

#   35 Tropical
#   32 Subtropical
#   22 Desert
#   13 Grassland
#   3 Temperate

koepClass<-c("Australia",22,13,32,3,35)
Biome<-c("Australia","Desert","Grassland","Subtropical","Temperate","Tropical")
kpColour<-c("black","darkred","darkgoldenrod1","darkmagenta","cornflowerblue","chartreuse4")

#Australia
tmin<-subset(final,Biome==koepClass[1])[c("tminMean","Month")]
tmax<-subset(final,Biome==koepClass[1])[c("tmaxMean","Month")]
prec<-subset(final,Biome==koepClass[1])[c("precMean","Month")]
plot(tmax$Month, tmax$tmaxMean,xlab="",ylab = "",ylim = c(10,40),
     main = Biome[1],type = "l",col="red")
mtext(substitute(paste(Temperature, B * degree, "C)"), list(B = " (")),
      side=2, line=2,cex=.8)
par(new = TRUE)
plot(prec$Month,prec$precMean,type = "l", axes = FALSE, bty = "n",
     xlab = "", ylab = "",col = "blue",ylim=c(5,110))
axis(side=4, at = pretty(5:100))
mtext("Prec", side=4, line=2,cex=.8)
kpdat<-subset(dat,select=DOY_PL)
par(new = TRUE)
hist(kpdat[,1],main="",xlab="",ylab = "",xlim = c(1,365),axes=FALSE,freq=FALSE,breaks=seq(1,366,length = 13))

for (k in 2:length(koepClass)){
  #plot tmin, tmax, and precip for biome
  tmin<-subset(final,Biome==koepClass[k])[c("tminMean","Month")]
  tmax<-subset(final,Biome==koepClass[k])[c("tmaxMean","Month")]
  prec<-subset(final,Biome==koepClass[k])[c("precMean","Month")]
  
  plot(tmax$Month, tmax$tmaxMean,xlab="",ylab = "",ylim = c(10,40),
       main = Biome[k],type = "l",col="red")
  mtext(substitute(paste(Temperature, B * degree, "C)"), list(B = " (")),
        side=2, line=2,cex=.8)
  mtext('Month', side=1, line=2,cex=.8)
  par(new = TRUE)
  plot(prec$Month,prec$precMean,type = "l", axes = FALSE, bty = "n",
       xlab = "", ylab = "",col = "blue",ylim=c(5,110))
  axis(side=4, at = pretty(5:100))
  mtext("Prec", side=4, line=2,cex=.8)
  kpdat<-subset(dat,koeppen==koepClass[k],select=DOY_PL)
  par(new = TRUE)
  hist(kpdat[,1],main="",xlab="",ylab = "",xlim = c(1,365),axes=FALSE,freq=FALSE,breaks=seq(1,366,length = 13))
} 

#legend("topright", legend=Biome, pch=1 , lwd=1,col=kpColour, bty="n", text.font=3)


#===============================================
#PDF of breeding density in 10 year increments since post 1970 and pre

#===============================================
dat<-read.csv(paste0(dat.dir,"PeakBreedigPointOfLayDayOfYear2015-07-12.csv"))

dat<-subset(dat,koeppen==3)
dat$Radians <-(dat$DOY_PL/366*360)*pi / 180

a<-subset(dat,year <1970)$DOY_PL
b<-subset(dat,year >=1970 & year< 1980)$DOY_PL
c<-subset(dat,year >=1980 & year< 1990)$DOY_PL
d<-subset(dat,year >=1990 & year< 2000)$DOY_PL
e<-subset(dat,year >=2000 & year< 2010)$DOY_PL
f<-subset(dat,year >=2010)$DOY_PL
plot(density(a),col="red")
lines(density(b),col="orange")
lines(density(c),col="yellow")
lines(density(d),col="green")
lines(density(e),col="blue")
lines(density(f),col="purple")



a<-subset(dat,year <1990)$Radians
b<-subset(dat,year >=1990)$Radians
c<-subset(dat,year >=1980 & year< 1990)$Radians
d<-subset(dat,year >=1990 & year< 2000)$Radians
e<-subset(dat,year >=2000 & year< 2010)$Radians
f<-subset(dat,year >=2010)$Radians




PLc<- circular(dat$Radians, units = "radians", zero=pi/2, 
               rotation = "clock")

#get midpoints of 12 months on circle
bins<-12
arc <- (2 * pi)/bins
sector <- seq(0, 2 * pi - (2 * pi)/bins, length = bins)
midpoints <- circular(seq(arc/2, 2 * pi - pi/bins, length = bins), units = "radians", zero=pi/2, 
                      rotation = "clock")

#plot(PLc, cex=1.1, bin=720, stack=TRUE, sep=0.035, shrink=1.8)
rose.diag(PLc, bins=12, col="darkgrey", cex=1.1, prop=1.7,axes = FALSE,ticks = FALSE, shrink=1.5)
#rose.diag(PLc, bins=12, col="darkgrey", cex=1.1, prop=1.3,rotation='clock',zero=pi/2,axes = FALSE,ticks = FALSE, shrink=1.8)
axis.circular(at=midpoints, labels=c(1:12), cex=1.1,tcl=0.075)
#ticks.circular(midpoints, rotation='clock', tcl=0.075)
#lines(density.circular(PLc, bw=40), lwd=2, lty=1)#Australia Density

kpColour<-c("chartreuse4","darkmagenta", "darkred","darkgoldenrod1","cornflowerblue")

date<-c("a","b","c","d","e","f") 

a<-circular(a, units = "radians", zero=pi/2, 
                rotation = "clock")
lines(density.circular(a, bw=40), lwd=3, lty=k,col = "red")
arrows.circular(mean.circular(a), lwd=3, lty=k,col = "red")

b<-circular(b, units = "radians", zero=pi/2, 
            rotation = "clock")
lines(density.circular(b, bw=40), lwd=3, lty=k,col = "orange")
arrows.circular(mean.circular(b), lwd=3, lty=k,col = "orange")

c<-circular(c, units = "radians", zero=pi/2, 
            rotation = "clock")
lines(density.circular(c, bw=40), lwd=3, lty=k,col = "yellow")
arrows.circular(mean.circular(c), lwd=3, lty=k,col = "yellow")

d<-circular(d, units = "radians", zero=pi/2, 
            rotation = "clock")
lines(density.circular(d, bw=40), lwd=3, lty=k,col = "green")
arrows.circular(mean.circular(d), lwd=3, lty=k,col = "green")

e<-circular(e, units = "radians", zero=pi/2, 
            rotation = "clock")
lines(density.circular(e, bw=40), lwd=3, lty=k,col = "blue")
arrows.circular(mean.circular(e), lwd=3, lty=k,col = "blue")

f<-circular(f, units = "radians", zero=pi/2, 
            rotation = "clock")
lines(density.circular(f, bw=40), lwd=3, lty=k,col = "purple")
arrows.circular(mean.circular(f), lwd=3, lty=k,col = "purple")


  
  


# =====================================
#clades of birds

#======================================

library(plyr)
library(Hmisc)
library(multcomp)
library(gplots)
#get data for clades
  dat.dir<-'/Users/daisy/Google Drive/PhD/Data/Observaitons/Cleaned/Breeding/'
  #read in observations
  finPL<-read.csv(paste0(dat.dir,'PeakBreedigPointOfLayDayOfYear2015-07-12.csv'))
  BP<-read.csv("/Users/daisy/Google Drive/PhD/BreedingTiming/tables/BreedingQuantiles2015-07-12.csv")
  #remove species not included in study
  #get traits
  inc<-read.csv('/Users/daisy/Google Drive/PhD/BreedingTiming/tables/SpeciesOfInterest_2015-07-15.csv')
  inc2<-subset(inc,remove==1)
  RMspecies<-inc2$ScientificName
  finPL<-finPL[finPL$Scientific.Name %nin% RMspecies,]
  clades<-as.data.frame(table(inc$Hackett.coarse.clades))
  keepCL<-subset(clades,Freq >= 5)$Var1
  badclades<-subset(clades,Freq < 5)
  inc3<-inc[inc$Hackett.coarse.clades %nin% badclades$Var1,]
  inc3<-subset(inc3,remove!=1)
  #subset bp to just Australia
  ausBP<-subset(BP, Region=="Australia")
  #merge data with inc
  inc3<-droplevels(merge(inc3,ausBP,by.x="ScientificName", by.y = "Species"))

#figure out what is going on with clades
  passer<-droplevels(subset(inc3, order == "Passeriformes"))
  nonPasser<-droplevels(subset(inc3, order != "Passeriformes"))
#at least 2 group means are different from each, 
#need to do more to find out which ones
#p < 0.05  
oneway.test(inc3$BreedingPeriod ~ inc3$Hackett.coarse.clades)
oneway.test(passer$BreedingPeriod ~ passer$Hackett.coarse.clades)
oneway.test(nonPasser$BreedingPeriod ~ nonPasser$Hackett.coarse.clades)

boxplot(passer$BreedingPeriod ~ passer$Hackett.coarse.clades)
boxplot(nonPasser$BreedingPeriod ~ nonPasser$Hackett.coarse.clades)
boxplot(inc3$BreedingPeriod ~ inc3$Hackett.coarse.clades)
boxplot(inc3$Quantile5 ~ inc3$Hackett.coarse.clades)
boxplot(inc3$Quantile95 ~ inc3$Hackett.coarse.clades)


# means and sd of breeding period by clades
ddply(inc3,~Hackett.coarse.clades,summarise,mean=mean(BreedingPeriod),sd=sd(BreedingPeriod))

#renames clades so they are A - J

inc3$clades2<-inc3$Hackett.coarse.clades
#  
# inc3$clades2<-rename(inc3$Hackett.coarse.clades,replace = c("Honeyeaters, Thornbills, Bristlebirds, Australasian Wrens, Allies"= "1",                                                                                                                              
#                                                   "Hawks, Eagles, Secretarybird" = "2", "Parrots" = "3",
#                                                   "Ducks, Geese, Other Waterfowl, Screamers"  = "4", "Ibises, Herons"  = "5",                                                                                                                                                                              
#                                                   "Cuckoo-Shrikes, Helmetshrikes, Shrike-Flycatchers, Vangas, Whipbirds, Shrike-Thrushes, Whistlers, Shrikes, Monarchs, Vireos, Orioles, Drongos, Fantails, Birds Of Paradise, Crows, Jays, Allies" = "6",
#                                                   "Waders, Allies" = "7","Gulls, Terns, Auks, Crab Plover" = "8",
#                                                   "Pigeons, Doves"= "9","Australasian Robins"="10"))


#abbreviated levels

levels(inc3$Hackett.coarse.clades) <- abbreviate(levels(inc3$Hackett.coarse.clades))

inc3$Hackett.coarse.clades <- with(inc3, reorder(Hackett.coarse.clades,
                                                 BreedingPeriod, mean))

lmClades<-lm(formula = BreedingPeriod ~ Hackett.coarse.clades, data = inc3)

anova(lmClades)

TukeyRegion2<- glht(lmClades, linfct=mcp(Hackett.coarse.clades="Tukey"))

lets <- cld(TukeyRegion2)$mcletters$Letters

library(doBy)
SE <- function(x)sd(x)/sqrt(length(x))
mn <- summaryBy(BreedingPeriod ~ Hackett.coarse.clades, data=inc3,
                FUN=c(mean,SE))



par(mar=c(10,4,4,4), cex.lab=0.7, las=2, cex.axis=0.7)
b <- with(mn, barplot2(BreedingPeriod.mean,
                  names.arg=Hackett.coarse.clades,
                  plot.ci=TRUE, ylim=c(0,330),
                  ci.l=BreedingPeriod.mean - BreedingPeriod.SE,
                  ci.u=BreedingPeriod.mean + BreedingPeriod.SE))
text(b, 300, lets)





summary(TukeyRegion2)

par(mar=c(510,4,1))
plot(TukeyRegion2)





plot(density(subset(passer,Hackett.coarse.clades==keepCL[1])$BreedingPeriod))

for (i in 2:length(keepCL)){
  
 lines(density(subset(inc3,Hackett.coarse.clades==keepCL[i])$BreedingPeriod))

}


####What I really want to know is which clades have shorter/longer breeding periods than other cladesddply(inc3,~Hackett.coarse.clades,summarise,mean=mean(BreedingPeriod),sd=sd(BreedingPeriod))




