rm(list = ls())
library(raster)
library(data.table)
library(dplyr)
#climate date
silo.dir<-"/Users/daisy/Google Drive/PhD/Data/Spatial/Climate/SILO/"
#function to round to 0.05
rndfunc<-function(x){
  return(round(x*2,1)/2)
}
#prepare locations
obs<-fread("/Users/daisy/Google Drive/PhD/Data/Observaitons/Cleaned/Breeding/PeakBreedingPointOfLayDayOfYear2015-10-17.csv",data.table=FALSE)
locs<-subset(obs, year>=1957,select=c("lat","lon"))
locs$lat<-rndfunc(locs$lat)
locs$lon<-rndfunc(locs$lon)
locs<-unique(locs)
#locations of monthly data
precip<-list.files(paste0(silo.dir,"Monthly_Precip/"),full.name=TRUE)
maxT<-list.files(paste0(silo.dir,"Monthly_Tmax/"),full.name=TRUE)
minT<-list.files(paste0(silo.dir,"Monthly_Tmin/"),full.name=TRUE)
#funciton to extract climate data
clim<-function(i){
  pre<-extract(raster(precip[i]),
               cbind(locs$lon,locs$lat))
  maxTemp<-extract(raster(maxT[i]),
                cbind(locs$lon,locs$lat))
  minTemp<-extract(raster(minT[i]),
                cbind(locs$lon,locs$lat))
  return(cbind(locs,pre,maxTemp,minTemp))
}
#empty list
obsClim<-list()
#go through each month from Jan 1957
for(i in 1:length(maxT)){
    climDF<-clim(i)
    obsClim[[i]]<-climDF
    message(i)
  }
#put in dataframe  
climDF2<-as.data.frame(do.call("rbind",obsClim)) 
#write to file  
write.csv(climDF2,paste0('/Users/daisy/Google Drive/PhD/BreedingTiming/tables/',
                         'SiloClimateLocationAllMonths', 
                         as.Date(Sys.Date()),'.csv'),row.names=FALSE)




