rm(list = ls())
library(raster)
library(data.table)
library(Hmisc)
library(lme4)
library(car)
library(multcomp)
library(reporttools)
library(doBy)
library(chemometrics)#***
library(boot)#***

obs.dir <-
  '/Users/daisy/Google Drive/PhD/Data/Observaitons/Cleaned/Breeding/'

# # Breeding quantiles
# BP<-read.csv("/Users/daisy/Google Drive/PhD/BreedingTiming/tables/BreedingQuantiles2015-10-17.csv")
# traits
inc <-
  fread(
    '/Users/daisy/Google Drive/PhD/BreedingTiming/tables/SpeciesOfInterest_2015-10-07.csv',data.table =
      FALSE
  )
#peak breeding period data
peak <-
  fread(paste0(obs.dir,"PeakBreedingPointOfLayDayOfYear2015-10-17.csv"),data.table =
          FALSE)
#make sure only have species interested in
inc2 <- subset(inc,remove == 1)
RMspecies <- inc2$ScientificName
peak <-
  peak[peak$Scientific.Name %nin% RMspecies,] #remove unwanted species

inc <-
  inc[c("ScientificName","Jetz.Patch.clade","Common.me","genus","family","order")]

#fix day and month
peak$day <-
  as.numeric(strftime(as.Date(peak$DOY_PL, origin = "1970-01-01"),format = "%d"))
peak$month <-
  as.numeric(strftime(as.Date(peak$DOY_PL, origin = "1970-01-01"),format = "%m"))
peak <- merge(peak,inc, by.x = "Scientific.Name",by.y = "ScientificName")
#peak$day<-formatC(peak$day,width=2,flag='0')
#peak$month<-formatC(peak$month,width=2,flag='0')

##reclassify koeppen zones
koeppen <-
  raster(
    paste0(
      '/Users/daisy/Google Drive/PhD/Data/Spatial/BOM_climate_zones/kpngrp_major/koepenReclassified.asc'
    ),
    proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84")
  )
#combine equatorial and tropical: 41 Equatorial, 35 Tropical
m <- c(35, 41, 35)
rclmat <- matrix(m, ncol = 3, byrow = TRUE)
koeppen <- reclassify(koeppen,rclmat)
#add koeppen zones to peak data
peak$koeppen <- extract(koeppen,data.frame(cbind(peak$lon, peak$lat)))
#convert data to radians to make circular vector
#peak$Radians <-(peak$DOY_PL/366*360)*pi / 180

#combine climate and peak obs
clim <-
  fread(
    paste0(
      "/Users/daisy/Google Drive/PhD/BreedingTiming/tables/",
      "SiloClimateLocationDate30DaySummary2015-11-11.csv"
    ),
    data.table = FALSE
  )
colnames(clim)[2] <- "meanTmax"
MClim <-
  fread(
    paste0(
      "/Users/daisy/Google Drive/PhD/BreedingTiming/tables/",
      "SiloClimateLocationAllMonths2015-11-03.csv"
    ),
    data.table = FALSE
  )

#add koeppen zones to MCLIM
locs <- as.data.frame(unique(cbind(MClim$lon,MClim$lat)))
locs$koeppen <- extract(koeppen,locs)
MClim2 <- merge(MClim,locs,by.x = c("lon","lat"),by.y = c("V1","V2"))
MClim2 <- na.omit(MClim2)
#remove white space
clim$day <-
  as.vector(apply(as.data.frame(clim[,"day"]),2,function(x)
    gsub('\\s+', '',x)))
clim$month <-
  as.vector(apply(as.data.frame(clim[,"month"]),2,function(x)
    gsub('\\s+', '',x)))
clim$ID <- with(clim,paste(ID,day,month,year))
peak$ID <- with(peak,(paste(lat,lon,day,month,year)))
peak <- subset(peak,year >= 1957 & year <= 2013)
climPK <- merge(peak,clim, by.x = "ID",by.y = "ID")

#Koeppen zones: 35-Tropical, 32-Subtropical, 22-Desert,13-Grassland,3-Temperate
koep <- c(35, 13, 32, 22,  3)

#make summary of relized vs. potential
Smry <- list()
for (i in 1:length(koep)) {
  #get means of observed
  obs <- subset(climPK,koeppen == koep[i])
  obsM <- as.matrix(obs[,c("sumPrec","meanTmin","meanTmax")])
  obsMean <- colMeans(obsM)
  #get means of potential
  ptnl <- subset(MClim2,koeppen == koep[i])
  ptnlM <- as.matrix(ptnl[,c("pre","minTemp","maxTemp")])
  ptnlMean <- colMeans(ptnlM)
  #run bootstrap for potential
  ptnlSummary <- list()
  for (ii in 1:10000) {
    a <- obsM[sample(1:nrow(obs),replace = TRUE),]
    ptnlSummary[[ii]] <- colMeans(a)
  }
  dfSum <- as.data.frame(do.call("rbind",ptnlSummary))
  #get data for table
  Smry[[i]] <- c(
    koep[i],
    colMeans(obsM),
    colMeans(ptnlM),
    quantile(dfSum$sumPrec, probs = c(0.025, 0.975)),
    quantile(dfSum$meanTmin, probs = c(0.025, 0.975)),
    quantile(dfSum$meanTmax, probs = c(0.025, 0.975))
  )
}


SMruDF<-do.call("rbind",Smry)


















