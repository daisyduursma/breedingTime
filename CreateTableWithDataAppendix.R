#Calculate the quantiles on circular data
#make graph of Australia wide data

rm(list = ls())
library(raster)
library(circular)
library(gtools)
library(reshape2)

#Species
# Breeding quantiles, use to get species
BP<-read.csv("/Users/daisy/Google Drive/PhD/BreedingTiming/tables/BreedingQuantiles2015-10-17.csv")

inc<-read.csv('/Users/daisy/Google Drive/PhD/BreedingTiming/tables/SpeciesOfInterest_2015-10-07.csv')
#keep only species of interest
inc<-subset(inc, remove !=1)
taxon<-inc[c("genus","family","ScientificName")]
BP<-merge(BP,taxon,all.x=TRUE,by.x="Species",by.y="ScientificName")
BP<-subset(BP, !is.na(BP$family))
#unique number of species
species<-as.character(unique(BP$Species))


CalcbinSize <- function(epochDates){
  d <- density(epochDates,na.rm = TRUE,from = 1, to = 366)
  ap <- approxfun(x=d$x, y=d$y)
  binSize<-integrate(ap, 1, 30.5)[[1]]
  for (mnth in 1:11){
    binSize<-c(binSize,integrate(ap, 30.5*mnth, (30.5*mnth)+30.5)[[1]])
  } 
  return(binSize)
}

# calcskew <- function(circulardates){
#   datesC<-circular( circulardates,modulo ="2pi", units="radians", rotation="counter")
#   Rbar<-rho.circular(datesC)#average clustering 0 is uncluseted 1 is all at same location
#   V<-1-Rbar#sample circular variance
#   t2t<- trigonometric.moment(datesC, p=2, center=TRUE)
#   bbar2 <- t2t$sin
#   skewness <- bbar2/(V**(3/2)) #skewness 
#   return(round(skewness,2))
# }
# 
#Koeppen zones
#   41 Equatorial
#   35 Tropical
#   32 Subtropical
#   22 Desert
#   13 Grassland
#   3 Temperate

#working directory
dat.dir<-'/Users/daisy/Google Drive/PhD/Data/Observaitons/Cleaned/Breeding/'
#read in observations
dat<-read.csv(paste0(dat.dir,'PointOfLayDayOfYear2015-10-07.csv'))
#figure out data types
dat<-droplevels(dat[dat$Scientific.Name %in% species,])


#egg
dat$type<- with(dat,ifelse((startEgg==startEgg & is.na(startHatch) & is.na(startYoung) & is.na(endBuild)),"solitary-egg",type))
#unknown
dat$type<- with(dat,ifelse((is.na(startEgg)&is.na(startEgg) & is.na(startHatch) & is.na(startYoung) & is.na(endBuild)),"unknown",type))
#young
dat$type<- with(dat,ifelse((is.na(startEgg)&is.na(startEgg) & is.na(startHatch) & !is.na(startYoung) & is.na(endBuild)),"young",type))
# multi
dat$type<- with(dat,ifelse((!is.na(startHatch)),"multi",type))     
dat$type<- with(dat,ifelse((!is.na(startEgg) &!is.na(startEgg) & !is.na(startYoung)),"multi",type))     
dat$type<- with(dat,ifelse((!is.na(startEgg) &!is.na(startEgg) & !is.na(endBuild)),"multi",type))
dat$type<- with(dat,ifelse((is.na(startEgg) &is.na(startEgg) & !is.na(endBuild)& !is.na(startYoung)),"multi",type))
#build plus unknown
dat$type<- with(dat,ifelse((is.na(startEgg) & is.na(startEgg) & !is.na(endBuild) & is.na(startYoung) & is.na(startHatch)& !is.na(startUnknown)),"unknown",type)) 


#add koeppen zones
koeppen<-raster(paste0('/Users/daisy/Google Drive/PhD/Data/Spatial/BOM_climate_zones/kpngrp_major/koepenReclassified.asc'),
                proj4string=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))
#combine equatorial and tropical: 41 Equatorial, 35 Tropical
m <- c(35, 41, 35)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
koeppen<- reclassify(koeppen, rclmat)


dat$koeppen<-extract(koeppen,data.frame(cbind(dat$lon, dat$lat)))
dat<-subset(dat,!is.na(koeppen))
#convert data to radians to make circular vector
dat$Radians <-(dat$DOY_PL/366*360)*pi / 180
#get list of incubation, fledging, clutch size etc
# inc<-read.csv('/Users/daisy/Google Drive/PhD/BreedingTiming/tables/SpeciesOfInterest_2015-10-07.csv')
# species<-inc$ScientificName
#species<-"Pomatostomus ruficeps"
#loop through species and convert to circular and find quantiles
#quantdat<-list()
peakBreeding <- list()
for (i in 1:length (species)) {
  spdat <- subset(dat,Scientific.Name == species[i])#get species data
  CommonName <-
    as.character(subset(inc,ScientificName == species[i],"Common.me",drop =
                          TRUE))
  Order <-
    as.character(subset(inc,ScientificName == species[i],order,drop = TRUE))
  #########################
  # koeppen zone analyses
  #########################
  #koeppen zones
  koep <- c(35,32,22,13,3)
  kpName <-
    c('Tropical','Subtropical','Desert','Grassland','Temperate')
  koepObs <- list()
  #colour<-c('darkgreen','green','chartreuse','darkgoldenrod4','purple','blue3')
  #koepQuant<-list()
  for (ii in 1:length(koep)) {
    kpdat <- as.numeric(subset(spdat,koeppen == koep[ii],Radians)[,1])
    if(length(kpdat)==0) {next}#if there are no observations go on to the next one
    kpobs <- subset(spdat,koeppen == koep[ii])
    # message(nrow(kpobs))
    if (length(kpdat) >= 50) {#more than 50 obs, calc peak
      #date of quantiles as long as over 50 obs in a region, number comes from Joys and Cricks
      kpquant <- quantile.circular(kpdat,c(.05,.5,.95),type = 8)
      kp5 <- round((kpquant[[1]] * 180 / pi) / 360 * 365)
      kp50 <- round((kpquant[[2]] * 180 / pi) / 360 * 365)
      kp95 <- round((kpquant[[3]] * 180 / pi) / 360 * 365)
      
      ###90% of data, remove obs outside of peak breeding
      kpObs90 <- if (kpquant[[1]] < kpquant[[3]]) {
        subset(spdat, Radians >= kpquant[[1]] &
                 Radians <= kpquant[[3]] & koeppen == koep[ii])
      }else
        rbind(
          subset(spdat, Radians >= kpquant[[1]] &
                   koeppen == koep[ii]),subset(spdat, Radians <= kpquant[[3]] &
                                                 koeppen == koep[ii])
        )
      kpObsExclude <- if (kpquant[[1]] < kpquant[[3]]) {
        rbind(
          subset(spdat, Radians < kpquant[[1]] &
                   koeppen == koep[ii]),subset(spdat, Radians > kpquant[[3]] &
                                                 koeppen == koep[ii])
        )
      }else
        subset(spdat, Radians < kpquant[[1]] &
                 Radians > kpquant[[3]] & koeppen == koep[ii])

      kpObs90$peak <- 1
      kpObsExclude$peak <- 0
      kpObsEnd <- smartbind(kpObs90,kpObsExclude)
     #message(nrow(kpObsEnd))
      
      koepObs[[ii]] <- kpObsEnd
      
    } else{
      if(as.numeric(nrow(kpobs))>0){
        kpobs$peak <- NA
        # message(nrow(kpobs))
        koepObs[[ii]] <-
        kpobs#if there was not 50 observations put the data back in
      }
    }
    rm(kpdat,kpobs,kpObs90,kpObsEnd,kpObsExclude,kpquant,kp5,kp50,kp95)
  }
  
  koepObs2 <- do.call("rbind",koepObs)
 # try(if(nrow(koepObs2)!=nrow(spdat)) stop("species counts not correct"))
  message(nrow(spdat))
  message(nrow(koepObs2))
  peakBreeding[[i]] <-
      do.call("rbind",koepObs)
  message(i)
}


# alldat<-na.omit(do.call("rbind",quantdat)[,c('Region','Species','CommonName','Order',
#                           'ObservationCount','Quantile5','Quantile50','Quantile95','BreedingPeriod','RbarAll','SkewAll','Rbar90','Skew90')])
peakBreedingAll<-do.call("rbind",peakBreeding)[,c("Scientific.Name","lat","lon","koeppen","DOY_PL","year","peak")]


############################

#add climate data


##############################

clim<-read.csv<-read.csv("/Users/daisy/Google Drive/PhD/BreedingTiming/tables/SiloClimateLocationDate30DaySummary2015-11-11.csv")
clim$month<-formatC(clim$month,width=2,flag='0')
clim$day<-formatC(clim$day,width=2,flag='0')

dat<-peakBreedingAll
dat$day<-as.numeric(strftime(as.Date(dat$DOY_PL, origin = "1970-01-01"),format = "%d"))
dat$month<-as.numeric(strftime(as.Date(dat$DOY_PL, origin = "1970-01-01"),format = "%m"))
dat$month<-formatC(dat$month,width=2,flag='0')
dat$day<-formatC(dat$day,width=2,flag='0')

dat$ID<-paste(dat$lat,dat$lon)
dat2<-merge(dat,clim,by=c("ID","day","month","year"),all.x=TRUE)#add climate data to table


#############################

#add ENSO Phase

#############################

SOI<-read.csv("/Users/daisy/Google Drive/PhD/Data/SOI/NOAA_OceanicNinoIndex5MonthPhase.csv")
colnames(SOI)<-c("Year",1:12)
SOI<-melt(SOI, id="Year")
SOI$month<-formatC(SOI$variable,width=2,flag='0') #format months 1 becomes 01
SOI$date<-paste0(SOI$Year,SOI$month) #make date string that matches NOAA SOI

dat2$date<-paste0(dat2$year,dat2$month)
dat3<-merge(dat2,SOI, by.x="date",by.y="date",all.x=TRUE)


##########################################
############ change names and keep column wanted ########
##########################################


dat4<-dat3[c("Scientific.Name","koeppen","DOY_PL","mean.maxT.","sumPrec","meanTmin","value","peak")]

colnames(dat4)<-c("Species","biome","DOY","Tmax","Prec","TMIN","SOI","peak")


write.csv(dat4,paste0('/Users/daisy/Google Drive/PhD/BreedingTiming/manuscript/TablesFigures/ManuscriptData',
          as.Date(Sys.Date()),'.csv'),row.names=FALSE)










