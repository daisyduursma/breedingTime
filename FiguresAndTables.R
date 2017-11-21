rm(list = ls())
library(raster)
library(circular)
#library(CircStats)
library(gtools)
library(Hmisc)

library(plyr)
library(Hmisc)
library(multcomp)
library(gplots)




obs.dir<-'/Users/daisy/Google Drive/PhD/Data/Observaitons/Cleaned/Breeding/'


#traits for species
#inc<-read.csv('/Users/daisy/Google Drive/PhD/BreedingTiming/tables/SpeciesOfInterest_2015-07-06.csv')
BP<-read.csv("/Users/daisy/Google Drive/PhD/BreedingTiming/tables/BreedingQuantiles2015-07-12.csv")
#   #remove species not included in study
#   #get traits
inc<-read.csv('/Users/daisy/Google Drive/PhD/BreedingTiming/tables/SpeciesOfInterest_2015-07-15.csv')
#
#get list of species not to include
inc2<-subset(inc,remove==1)
RMspecies<-inc2$ScientificName

#peak breeding period data
peak<-read.csv(paste0(obs.dir,"PeakBreedigPointOfLayDayOfYear2015-07-14.csv"))
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


# =====================================

# Figure 1. circular plot of Australia with density curves for biomes, 2X1 graph with UK on other side

# ===================================== 
{
  pdf(file = "/Users/daisy/Google Drive/PhD/BreedingTiming/figures/breedingObservationDensity20150929.pdf",
       width = 6, height = 3)
  
par(mfrow = c(1,3),
      oma = c(4,2,0,0),xpd=NA)
par(mar = c(1,0,0,0))    
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
axis.circular(at=midpoints, labels=c(1:12), cex=1.1,tcl=0.075)
#Koeppen zones: 35-Tropical, 32-Subtropical, 22-Desert,3-Grassland,3-Temperate
koepClass<-c(35,32,22,13,3)
kpColour<-c("green4","darkseagreen", "coral1","darkgoldenrod1","cornflowerblue")
#plot kernal density lines for each koeppen region
lines<-c(1,2,4,5,6)
for (k in 1:length(koepClass)){
  kpdat<-circular(subset(peak,koeppen==koepClass[k],select = Radians), units = "radians", zero=pi/2, 
                  rotation = "clock")

  lines(density.circular(kpdat, bw=12), lwd=1, lty=lines[k],
        col = kpColour[k],shrink=.6)
  #arrows.circular(median(kpdat), lwd=2, lty=lines[k],col = kpColour[k])#arrow pointing to median date
}

#=======AUS Number of species breeding in any month
#get unique species per month in Australia and make rose diagram
monthdat<-peak[!duplicated(peak[c("Scientific.Name","month")]),]
Mthc<- circular(monthdat$Radians, units = "radians", zero=pi/2, 
               rotation = "clock")

rose.diag(Mthc, bins=12, col="grey90", cex=1.1, prop=1.7,
          axes = FALSE,ticks = FALSE, shrink=1.5)
axis.circular(at=midpoints, labels=c(1:12), cex=1.1,tcl=0.075)

#for biomes
koepClass<-c(35,32,22,13,3)
for (k in 1:length(koepClass)){
  kpdat<-subset(monthdat,koeppen==koepClass[k])
  kpdat<- kpdat[!duplicated(kpdat[c("Scientific.Name","month")]),]
  kpdat<-circular(kpdat$Radians, units = "radians", zero=pi/2, 
                  rotation = "clock")

  lines(density.circular(kpdat, bw=12), lwd=1,
        lty=lines[k],col = kpColour[k],shrink=.6)
  #arrows.circular(median(kpdat), lwd=2, lty=lines[k],col = kpColour[k])
}
#=======UK Number of species breeding in any month
#par(mar = c(1,0,0,0)) 
#get UK data
UK<-read.csv('/Users/daisy/Google Drive/PhD/BreedingTiming/tables/Joys_etal_Breeding_periods_England.csv')

#for each species calculate the months breeding
 
spDay<-list()
for(sp in 1:length(UK$Species)){
  spdat<-UK[sp,]
  per5<-with(spdat,ifelse(is.na(PercentilesAllYears.5),PercentilesPost1990.5,PercentilesAllYears.5))
  per95<-with(spdat,ifelse(is.na(PercentilesAllYears.95),PercentilesPost1990.95,PercentilesAllYears.95))
  spdat<-data.frame(seq(per5,per95),paste(spdat$Species,spdat$species.2))
  spdat$Radians <-(spdat[,1]/366*360)*pi / 180
  spDay[[sp]]<-spdat
}
UKDay<-do.call("rbind",spDay)
colnames(UKDay)<-c("day","species","Radians")

CalcEpochMonth <- function(PointLay){
  mn<-as.numeric(strftime(as.Date(PointLay, origin = "1970-01-01"),format = "%m"))
  return(mn)
}
UKDay$month<-CalcEpochMonth(UKDay$day)
UKmonthdat<-UKDay[!duplicated(UKDay[c("species","month")]),]

Mthc<- circular(UKmonthdat$Radians, units = "radians", zero=pi/2, 
                rotation = "clock")


rose.diag(Mthc, bins=12, col="grey90", cex=1.1, prop=1.7,
          axes = FALSE,ticks = FALSE, shrink=1.5)

axis.circular(at=midpoints, labels=c(1:12), cex=1.1,tcl=0.075)

#lines(density.circular(Mthc, bw=10), lwd=3,col = "black")

kpColour<-c("green4","darkseagreen", "coral1","darkgoldenrod1","cornflowerblue")

legend(c("Desert","Grassland"),
         pch=1 , lwd=1, bty="n", text.font=1, 
         col = c("coral1","darkgoldenrod1"),x = -8, y =-1.35,cex=1.1)
legend(c("Subtropical","Temperate"),
       pch=1 , lwd=1, bty="n", text.font=1, 
       col = c("darkseagreen","cornflowerblue"),x = -4.5, y = -1.35,cex=1.1)
legend(c("Tropical"),
       pch=1 , lwd=1, bty="n", text.font=1, 
       col = c("green4"),x = -1, y = -1.35,cex=1.1)
 text("a",x=-8,y=1.15,font=2,cex=1.1)
 text("b",x=-4.6,y=1.15, font =2,cex=1.1)
text("c",x=-1.35,y=1.15, font =2,cex=1.1)

dev.off()
}

##################################

#PDFs of observation years
##################################
pdf(file = "/Users/daisy/Google Drive/PhD/BreedingTiming/figures/ObservationYear20150907.pdf",
    width = 6, height = 4)


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
plot(density(altYear), main ="PDF year of collection",xlab="Year",xlim=c(1770,2015),
     col="darkorchid4", lwd=2,cex.lab=1)
lines(density(NRSYear), col="azure4", lwd=2)
lines(density(musYear), col="darkblue", lwd=2)
lines(density(ABBBSyear), col="darkred", lwd=2)
lines(density(ebirdyear), col="darkgoldenrod2", lwd=2)
legend("topleft",c("Atlas","NRS","Museum","ABBBS","eBIRD"),text.col=c("darkorchid4","azure4","darkblue","darkred","darkgoldenrod2"),cex=1)

dev.off()


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

# peak<-read.csv(paste0(dat.dir,"PeakBreedigPointOfLayDayOfYear2015-07-12.csv"))
# peak<-peak[peak$Scientific.Name %nin% RMspecies,]
nrow(peak)



#==========================

# Figure 2. temperature and observations by biome

#==========================

library(gtools)

#climate variabls
clim<-c('prec','tmax','tmin' )
climateDir<-"/Users/daisy/Google Drive/PhD/Data/Spatial/Climate/EMASTClimate_mmn"
fPARDir<-"/Users/daisy/Google Drive/PhD/Data/Spatial/fPAR/"
#koeppen zones, combine equatorial and tropic
koeppen<-raster(paste0('/Users/daisy/Google Drive/PhD/Data/Spatial/BOM_climate_zones/kpngrp_major/koepenReclassified.asc'),
                proj4string=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))
m <- c(35, 41, 35)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
koeppen<- reclassify(koeppen, rclmat)

r1 <- raster("/Users/daisy/Google Drive/PhD/Data/Spatial/Climate/EMASTClimate_mmn/tmin/eMAST_ANUClimate_mmn_tmin_v1m0_01.nc")
r2<-raster("/Users/daisy/Google Drive/PhD/Data/Spatial/fPAR/Aust_ftot_8km_monthly_1980-2011_v5_AVG_month12.asc")
#make koeppen zones (more course than climate) so that they have same
#resolution and extent as climate
koeppen<-resample(koeppen,r1,method="ngb")
kpfpar<-resample(koeppen,r2,method="ngb")

#Peak breeding observations
Locs <- SpatialPoints(cbind(peak$lon, peak$lat),
                      proj4string = CRS('+proj=longlat +datum=WGS84'))

#climatic conditions for biomes
endClimDat<-list()

for (c in 1:length(clim)) {
  #get climate files
  files <- list.files(paste0(climateDir,"/",clim[c]))
  # make a loop to sort through the files
  climdat<-list()
  for (i in 1:length(files)) {
    # make a raster for the desired time period
    r1 <- raster(paste0(climateDir,"/",clim[c],"/",files[i]))
    mr1<-cellStats(r1,mean)
    sdr1<-cellStats(r1,sd)
    month<-strsplit(strsplit(files[i],"_")[[1]][6],".nc")[[1]][1]
    # extract the data and write it directly to the csv file
    MeanBiomes <- zonal(r1,koeppen, fun = mean,na.rm=TRUE)
    sdBiomes <- zonal(r1,koeppen, fun = sd,na.rm=TRUE)
    
    #fPAR
    fPAR<-raster(paste0(fPARDir,"Aust_ftot_8km_monthly_1980-2011_v5_AVG_month",month,".asc"))
    MeanfPAR <- zonal(fPAR,kpfpar, fun = mean,na.rm=TRUE)
    sdfPAR<-zonal(fPAR,kpfpar, fun = sd,na.rm=TRUE)
    mfPARAUS<-cellStats(raster::mask(fPAR,kpfpar),mean)
    sdfPARAUS<-cellStats(raster::mask(fPAR,kpfpar),sd)
    
    # reasign the column name to the one we generated
    enddat<-merge(MeanBiomes,sdBiomes, by="zone")
    enddat<-merge(enddat,MeanfPAR, by="zone")
    enddat<-merge(enddat,sdfPAR, by="zone")
    colnames(enddat)<-c("Biome",paste0(clim[c],"Mean"),paste0(clim[c],"Sd"),"MeanfPAR","sdfPAR")
    enddat<-rbind(enddat,c("Australia",mr1,sdr1, mfPARAUS,sdfPARAUS))
    enddat$Month<-month
    climdat[[i]]<-enddat
   
  }
  endClimDat[[c]]<-data.frame(do.call("smartbind",climdat))
}
#make table of summary  
final<-do.call("cbind",endClimDat)[c("Biome","Month","precMean","precSd","tmaxMean","tmaxSd",
                                     "tminMean","tminSd","MeanfPAR","sdfPAR")]


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
fp<-subset(final,Biome==koepClass[1])[c("MeanfPAR","Month")]
plot(tmax$Month, tmax$tmaxMean,xlab="",ylab = "",ylim = c(10,40),
     main = Biome[1],type = "l",col="red")
mtext(substitute(paste(Temperature, B * degree, "C)"), list(B = " (")),
      side=2, line=2,cex=.8,adj=0)
par(new = TRUE)
plot(prec$Month,prec$precMean,type = "l", axes = FALSE, bty = "n",
     xlab = "", ylab = "",col = "blue",ylim=c(5,110))
axis(side=4, at = pretty(5:100))
mtext("Prec", side=4, line=2,cex=.8)
par(new = TRUE)
plot(fp$Month,fp$MeanfPAR,type = "l", axes = FALSE, bty = "n",
     xlab = "", ylab = "",col = "brown",ylim=c(0,.7),lty = 3,lwd=2)
kpdat<-subset(peak,select=DOY_PL)
par(new = TRUE)
hist(kpdat[,1],main="",xlab="",ylab = "",xlim = c(1,365),axes=FALSE,freq=FALSE,breaks=seq(1,366,length = 13))

for (k in 2:length(koepClass)){
  #plot tmin, tmax, and precip for biome
  tmin<-subset(final,Biome==koepClass[k])[c("tminMean","Month")]
  tmax<-subset(final,Biome==koepClass[k])[c("tmaxMean","Month")]
  prec<-subset(final,Biome==koepClass[k])[c("precMean","Month")]
  fp<-subset(final,Biome==koepClass[k])[c("MeanfPAR","Month")]
  
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
  kpdat<-subset(peak,koeppen==koepClass[k],select=DOY_PL)
  par(new = TRUE)
  hist(kpdat[,1],main="",xlab="",ylab = "",xlim = c(1,365),axes=FALSE,freq=FALSE,breaks=seq(1,366,length = 13))
  par(new = TRUE)
  plot(fp$Month,fp$MeanfPAR,type = "l", axes = FALSE, bty = "n",
       xlab = "", ylab = "",col = "brown",ylim=c(0,.7),lty = 3,lwd=2)
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
  #dat.dir<-'/Users/daisy/Google Drive/PhD/Data/Observaitons/Cleaned/Breeding/'
  #read in observations
  #finPL<-read.csv(paste0(dat.dir,'PeakBreedigPointOfLayDayOfYear2015-07-12.csv'))
   inc2<-subset(inc,remove==1)
#   RMspecies<-inc2$ScientificName
#   finPL<-finPL[finPL$Scientific.Name %nin% RMspecies,]
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
#  passer<-droplevels(subset(inc3, order == "Passeriformes"))
#  nonPasser<-droplevels(subset(inc3, order != "Passeriformes"))
#at least 2 group means are different from each, 
#need to do more to find out which ones
#p < 0.05  
# oneway.test(inc3$BreedingPeriod ~ inc3$Hackett.coarse.clades)
# oneway.test(passer$BreedingPeriod ~ passer$Hackett.coarse.clades)
# oneway.test(nonPasser$BreedingPeriod ~ nonPasser$Hackett.coarse.clades)
# 
# boxplot(passer$BreedingPeriod ~ passer$Hackett.coarse.clades)
# boxplot(nonPasser$BreedingPeriod ~ nonPasser$Hackett.coarse.clades)
# boxplot(inc3$BreedingPeriod ~ inc3$Hackett.coarse.clades)
# boxplot(inc3$Quantile5 ~ inc3$Hackett.coarse.clades)
# boxplot(inc3$Quantile95 ~ inc3$Hackett.coarse.clades)


# means and sd of breeding period by clades
summary<-ddply(inc3,~Hackett.coarse.clades,summarise,mean=mean(BreedingPeriod),sd=sd(BreedingPeriod))

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



# 
# 
# 
# 
# 
# 
# 
# summary(TukeyRegion2)
# 
# par(mar=c(510,4,1))
# plot(TukeyRegion2)
# 
# 
# plot(density(subset(passer,Hackett.coarse.clades==keepCL[1])$BreedingPeriod))
# 
# for (i in 2:length(keepCL)){
#   
#  lines(density(subset(inc3,Hackett.coarse.clades==keepCL[i])$BreedingPeriod))
# 
# }


####What I really want to know is which clades have shorter/longer breeding periods than other cladesddply(inc3,~Hackett.coarse.clades,summarise,mean=mean(BreedingPeriod),sd=sd(BreedingPeriod))

#=====================

#Is there differences between thermal niche of biomes?

#==================================================

#zf<-subset(peak, !is.na(koeppen)&Scientific.Name=="Taeniopygia guttata")
koepPeak<-subset(peak, !is.na(koeppen))

#means of species in biomes
koepPeak<-ddply(koepPeak, .(Scientific.Name,koeppen), summarize,  mmtmax=mean(mmtmax))
#means of biomes weigthed by biome
b<-ddply(a, .(koeppen), summarize,  mmtmax=mean(mmtmax))



levels(koepPeak$koeppen) <- abbreviate(levels(koepPeak$koeppen))
koepPeak$koeppen <- with(koepPeak, reorder(koeppen,
                                   mmtmax, mean))
lmClades<-lm(formula = mmtmax ~ koeppen, data = koepPeak)
anova(lmClades)
TukeyRegion2<- glht(lmClades, linfct=mcp(koeppen="Tukey"))
lets <- cld(TukeyRegion2)$mcletters$Letters
library(doBy)
SE <- function(x)sd(x)/sqrt(length(x))
mn <- summaryBy(mmtmax ~ koeppen, data=koepPeak,
                FUN=c(mean,SE))
par(mar=c(10,4,4,4), cex.lab=0.7, las=2, cex.axis=0.7)
b <- with(mn, barplot2(mmtmax.mean,
                       names.arg=koeppen,
                       plot.ci=TRUE, ylim=c(0,40),
                       ci.l=mmtmax.mean - mmtmax.SE,
                       ci.u=mmtmax.mean + mmtmax.SE))
text(b, 35, lets)


library(plyr)
#mean of observatons


#mean of species
#=====================================
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

EGM_ALL<-list()
for ( i in 1: length(unique(peak$Scientific.Name))){
  sp<-subset(peak,Scientific.Name==unique(peak$Scientific.Name)[[i]]) 
  obs<-as.data.frame(table(sp$month))$Freq
  EGM<-as.data.frame(CalcEGM(obs))
  colnames(EGM)<-"AUS"
  
  for(ii in 1:length(unique(peak$koeppen))){
    spbiome<-subset(sp,koeppen==unique(peak$koeppen)[[ii]])
    obsBiome<-as.data.frame(table(spbiome$month))$Freq
    obsBiome<- if (length(obsBiome)==12) { obsBiome
    } else c(obsBiome,rep(0,12-length(obsBiome)))
    EGM$biome<-ifelse(nrow(spbiome)>=50,CalcEGM(obsBiome),NA)
    colnames(EGM)[ncol(EGM)]<-paste0("biome_",unique(peak$koeppen)[[ii]])
  }
  
  EGM$species<-unique(peak$Scientific.Name)[[i]]
  
  EGM_ALL[[i]]<-EGM
}
finEGM<-as.data.frame(do.call("rbind",EGM_ALL))

finEGM<-finEGM[,c("species","AUS","biome_13","biome_3","biome_32","biome_35","biome_22")]

#Koeppen zones

#   35 Tropical and Equatorial
#   32 Subtropical
#   22 Desert
#   13 Grassland
#   3 Temperate

colMeans(finEGM[,2:7],na.rm=TRUE)


require(reshape2)
finEGM2<-melt(finEGM,id = "species")
finEGM2<-na.omit(finEGM2)

finEGM2$variable <- with(finEGM2, reorder(variable,
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

#==================================

#see if anything is going on in south west Australia, 
#where there are warmer winter and spring (minimum tempertures) and 25% drier

#======================================
#get data
#koepPeak<-subset(peak, !is.na(koeppen))

#pick out species
#keeps obs south of -26 and west of 125 with at least 100 obs
SWobs<-subset(peak, lat <= -26 & lon <=125)
spSum<-as.data.frame(table(SWobs$Scientific.Name))
spSum<-as.vector(droplevels(subset(spSum,Freq>=200)$Var1))
SWobs<-SWobs[SWobs$Scientific.Name %in% spSum,] #
ausBP<-subset(BP, Region=="Australia")
inc_new<-droplevels(merge(inc,ausBP,by.x="ScientificName", by.y = "Species"))

#make figures of significant differences and new dataframe with modified point of lay
newSWobs<-list()
for (i in 1:length(spSum)){
  start<-subset(inc_new,ScientificName==spSum[i],Quantile5)
   d1<-subset(SWobs,Scientific.Name==spSum[i] &year<1990)
   d1$type<-"pre"
   d2<-subset(SWobs,Scientific.Name==paste(spSum[i])&year>=1990)
   d2$type<-"post"
   d3<-rbind(d1,d2)
   #adjust dates so that everything less than start is added to 365, then do not need to worry about 1 being equally as close as 364
   d3$MODIFIED_DOY_PL<-with(d3,ifelse(DOY_PL<start$Quantile5,DOY_PL+365,DOY_PL))
   levels(d3$type) <- abbreviate(levels(d3$type))
    d3$type <- with(d3, reorder(type,
                              MODIFIED_DOY_PL, mean))
   lmClades<-lm(formula = MODIFIED_DOY_PL ~ type, data = d3)
   anova(lmClades)
   TukeyRegion2<- glht(lmClades, linfct=mcp(type="Tukey"))
   lets <- cld(TukeyRegion2)$mcletters$Letters
   library(doBy)
   SE <- function(x)sd(x)/sqrt(length(x))
   mn <- summaryBy(MODIFIED_DOY_PL ~ type, data=d3,
                   FUN=c(mean,SE))
   par(mar=c(10,4,4,4), cex.lab=0.7, las=2, cex.axis=0.7)
   b <- with(mn, barplot2(MODIFIED_DOY_PL.mean,
                          names.arg=type,
                          main=spSum[i],
                          plot.ci=TRUE, ylim=c(0,350),
                          ci.l=MODIFIED_DOY_PL.mean - MODIFIED_DOY_PL.SE,
                          ci.u=MODIFIED_DOY_PL.mean + MODIFIED_DOY_PL.SE))
   text(b, 150, lets)
   text(b,100,paste0("pre1990 = ", nrow(d1), ", post1990 = ", nrow(d2)))
   d1NRS<-subset(d1,sourceName == "NestRecordScheme")
   d2NRS<-subset(d2,sourceName == "NestRecordScheme")
   text(b,50,paste0("NRS: pre1990 = ", nrow(d1NRS), ", post1990 = ", nrow(d2NRS)))
   newSWobs[[i]]<-d3
}
newSWobs2<-do.call("rbind",newSWobs)


#species of interest with at least 200 obs and more than 50 in pre and post, show temperal change
SI<-c("Acanthiza chrysorrhoa", 
      "Anthochaera carunculata",
      "Corvus coronoides",
      "Cracticus tibicen",
      "Cygnus atratus",
      "Eolophus roseicapillus",
      "Petrochelidon ariel",
      "Phylidonyris novaehollandiae",
      "Rhipidura fuliginosa",
      "Rhipidura leucophrys",
      "Zosterops lateralis")

#keep only species with significant change
SWobs2<-newSWobs2[newSWobs2$Scientific.Name %in% SI,]
#divide the data into two groups and make tukey figure
d1<-subset(SWobs2,year<1990)
d1$dates<-"pre"
d2<-subset(SWobs2,year>=1990)
d2$dates<-"post"
SWobs2<-rbind(d1,d2)

lmSW<-lm(formula = MODIFIED_DOY_PL ~ dates*Scientific.Name, data = SWobs2)
anova(lmSW)
with(SWobs2, interaction.plot(dates,Scientific.Name,MODIFIED_DOY_PL))
TukeyRegion2<- glht(lmClades, linfct=mcp(type="Tukey"))
lets <- cld(TukeyRegion2)$mcletters$Letters
library(doBy)
SE <- function(x)sd(x)/sqrt(length(x))
mn <- summaryBy(DOY_PL ~ type, data=d3,
                FUN=c(mean,SE))
par(mar=c(10,4,4,4), cex.lab=0.7, las=2, cex.axis=0.7)
b <- with(mn, barplot2(DOY_PL.mean,
                       names.arg=type,
                       main=spSum[i],
                       plot.ci=TRUE, ylim=c(0,300),
                       ci.l=DOY_PL.mean - DOY_PL.SE,
                       ci.u=DOY_PL.mean + DOY_PL.SE))
text(b, 150, lets)

##############################
#make table of 5 year means and standard errors, using this to make figures like in Crick's Nature paper
##############################
SWsummary<-list()
for (si in 1:length(SI)){
  #functions to calculate standard error of mean
  SEc <- function(x)round((sd.circular(x)/sqrt(length(x))*180/pi)/360*365)
  MNc <- function(x)round((quantile.circular(x,c(.5),type=8)[[1]]*180/pi)/360*365)
  
  #data cut into chunks
  obsPre65<-subset(SWobs2,Scientific.Name==SI[si]&year<1965,Radians)
  obs70<-subset(SWobs2,Scientific.Name==SI[si]&year<1975 & year>=1965,Radians)
  obs80<-subset(SWobs2,Scientific.Name==SI[si]&year<1985 & year>=1975,Radians)
  obs90<-subset(SWobs2,Scientific.Name==SI[si]&year<1995 & year>=1985,Radians)
  obs00<-subset(SWobs2,Scientific.Name==SI[si]&year<2005 & year>=1995,Radians)
  obs10<-subset(SWobs2,Scientific.Name==SI[si]&year<2015 & year>=2005,Radians)
 #means and standard errors
  summary<-as.data.frame(rep(SI[si],6))
  colnames(summary)<-"species"
  summary$mean<-c(
             MNc(obsPre65),
             MNc(obs70),
             MNc(obs80),
             MNc(obs90),
             MNc(obs00),
             MNc(obs10))
  summary$SE<-c(SEc(obsPre65),
             SEc(obs70),
             SEc(obs80),
             SEc(obs90),
             SEc(obs00),
             SEc(obs10))
  summary$year<-c(1965,1970,1980,1990,2000,2010)
  SWsummary[[si]]<-summary
  }

a<-do.call("rbind",SWsummary)

plot(mean~year,data=b,ylim=c(150,320),type="l", col="blue",ylab="Mean lay date of first eggs")

for (si in 1:length(SI)-1){
  b<-subset(a, species ==SI[si])

points(mean~year,data=b,type="l", col="blue")
  }

plot(mean~year,data=b,ylim=c(200,365),type="l", col="blue",ylab="Mean lay date of first eggs")
plot(mean~year,data=b,ylim=c(200,365),type="l", col="blue",ylab="Mean lay date of first eggs")
plot(mean~year,data=b,ylim=c(200,365),type="l", col="blue",ylab="Mean lay date of first eggs")



###################################
#SOI figure

##################################

#temperate region Breeding Period
#Koeppen zones: 35-Tropical, 32-Subtropical, 22-Desert,3-Grassland,3-Temperate
peakTEMP<-subset(peak, koeppen ==3)
peakTEMP$month<-formatC(peakTEMP$month,width=2,flag='0')
peakTEMP$date<-paste0(peakTEMP$year,peakTEMP$month)
peakTEMP<-subset(peakTEMP, date >=195101)

#read in Southern Oscillation Index
SOI<-read.csv("/Users/daisy/Google Drive/PhD/Data/SOI/NOAA_SouthernOscillationIndex.csv")
peakTEMP<-merge(peakTEMP,SOI, by.x="date",by.y="Date",all.x=TRUE)
#get el Nino obs and species with more than 100 observations
eN<-subset(peakTEMP,Value>=.5)
spSumeN<-as.data.frame(table(eN$Scientific.Name))
sp_eN<-as.vector(droplevels(subset(spSumeN,Freq>=200)$Var1))
eN<-eN[eN$Scientific.Name %in% spSum_eN,] #
#get la Nina obs and species with more than 100 obs
lN<-subset(peakTEMP,Value<=-.5)
spSumlN<-as.data.frame(table(lN$Scientific.Name))
sp_lN<-as.vector(droplevels(subset(spSumlN,Freq>=200)$Var1))
lN<-lN[lN$Scientific.Name %in% spSum_lN,] #

#name of species in both, limit data to these species
fin_sp<-intersect(sp_eN,sp_lN)
eN<-eN[eN$Scientific.Name %in% fin_sp,] #
eN<-eN[eN$Scientific.Name %in% spSum_eN,] #

hist(lN$DOY_PL,main="",xlab="",ylab = "",xlim = c(1,365),axes=FALSE,freq=FALSE,breaks=seq(1,366,length = 13))
hist(eN$DOY_PL,main="",xlab="",ylab = "",xlim = c(1,365),axes=FALSE,freq=FALSE,breaks=seq(1,366,length = 13))











##############################

#temperate vs Tropic

##############################

# three species having big differneces
SW<-c("Pardalotus striatus","Charadrius ruficapillus","Taeniopygia bichenovii")

#temperate region Breeding Period
tempBP<-subset(BP, Region =="Temperate")
#Tropic region Breeding Period
tropBP<-subset(BP, Region =="Tropical")
both<-merge(tempBP,tropBP,by = "Species")
both$mediandiff<- with(both,Quantile50.x-Quantile50.y)
both$startdiff<-with(both,Quantile5.x-Quantile5.y)
both$enddiff<-with(both,Quantile95.x-Quantile95.y)

SWdat<-both[both$Species %in% SW,]


#figure of 3 birds breeding

#pdf(file = "/Users/daisy/Google Drive/PhD/BreedingTiming/figures/breedingObservationDensity20150907.pdf",
 #   width = 6, height = 3)

par(mfrow = c(1,3),
    oma = c(4,0,0,0) + 0.1,xpd=NA)
par(mar = c(1,3,0,0))    
#get midpoints of 12 months on circle
bins<-12
arc <- (2 * pi)/bins
sector <- seq(0, 2 * pi - (2 * pi)/bins, length = bins)
midpoints <- circular(seq(arc/2, 2 * pi - pi/bins, length = bins), units = "radians", zero=pi/2, 
                      rotation = "clock")
#Koeppen zones: 35-Tropical, 3-Temperate
koepClass<-c(35,3)
kpColour<-c("green4","cornflowerblue")

#plot kernal density lines for each koeppen region
lines<-c(1,2,4,5,6)

for (SWsp in 1:length(SW)){
  rose.diag(PLc, bins=12, col="grey90", cex=1.1, prop=1.7,
            axes = FALSE,ticks = FALSE, shrink=1.2)

  #add month labels
  axis.circular(at=midpoints, labels=c(1:12), cex=1.1,tcl=0.075)
  
for (k in 1:length(koepClass)){
  

  kpdat<-circular(subset(peak,koeppen==koepClass[k] & Scientific.Name == SW[SWsp],
                  select = Radians), 
                  units = "radians", 
                  zero=pi/2, 
                  rotation = "clock")
  lines(density.circular(kpdat, bw=12), lwd=2, lty=lines[k],col = kpColour[k])

  mtext(SW[SWsp], side=2, line=1, at=1, cex=.8)
  
}

}









