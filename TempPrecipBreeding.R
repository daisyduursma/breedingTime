rm(list = ls())
library(raster)
library(data.table)
library(dplyr)

silo.dir<-"/Users/daisy/Google Drive/PhD/Data/Spatial/Climate/SILO/"

#climate data already created

climdat<-fread(paste0("/Users/daisy/Google Drive/PhD/BreedingTiming/tables/",
                   "SiloClimateLocationDate30DaySummary2015-10-26.csv"),
            data.table = FALSE)

climdat$day<-as.vector(apply(as.data.frame(climdat[,"day"]),2,function(x)gsub('\\s+', '',x)))
climdat$month<-as.vector(apply(as.data.frame(climdat[,"month"]),2,function(x)gsub('\\s+', '',x)))



precip<-list.files(paste0(silo.dir,"Daily_Precip/"))
maxT<-list.files(paste0(silo.dir,"Daily_Tmax/"))
minT<-list.files(paste0(silo.dir,"Daily_Tmin/"))

obs<-fread("/Users/daisy/Google Drive/PhD/Data/Observaitons/Cleaned/Breeding/PeakBreedingPointOfLayDayOfYear2015-10-17.csv",data.table=FALSE)
#calculate day of month
obs$day<-as.numeric(strftime(as.Date(obs$DOY_PL, origin = "1970-01-01"),format = "%d"))
obs$month<-as.numeric(strftime(as.Date(obs$DOY_PL, origin = "1970-01-01"),format = "%m"))
#keep everything after 1970
obs<-subset(obs,year>=1957&year<=2013)
obs<-obs[,c("lat","lon","year","month","day")]
obs<-obs[!duplicated(obs),]
obs2013<-subset(obs,year==2013&month<=10)
obs<-subset(obs,year<2013)
obs<-rbind(obs,obs2013)
obs1957m1<-subset(obs,year==1957&month==1&day>=15)
obs1957other<-subset(obs,year==1957&month>1)
obs<-subset(obs,year>1957)
obs<-rbind(obs,obs1957m1,obs1957other)



#function for returning prec, and tmax
clim<-function(day,lon,lat){
  dayD<-strftime(day,format = "%d")
  dayM<-strftime(day,format = "%m")
  dayY<-strftime(day,format = "%Y")
  pre<-extract(raster(paste0(silo.dir,"Daily_Precip/rai",dayY,dayM,dayD,".tif")),
               cbind(lon,lat))
  maxT<-extract(raster(paste0(silo.dir,"Daily_Tmax/Tmax",dayY,dayM,dayD,".tif")),
                cbind(lon,lat))
  minT<-extract(raster(paste0(silo.dir,"Daily_Tmin/Tmin",dayY,dayM,dayD,".tif")),
                cbind(lon,lat))
  return(cbind(lat,lon,pre,maxT,minT))
}



#get unique days
dayU<-with(obs,unique(paste(year,month,day)))
obsdone<-with(climdat,unique(paste(year,month,day)))
dayU<-setdiff(dayU,obsdone)



obsClim<-list()
for(i in 2975:length(dayU)){
  date<-strsplit(dayU[i]," ")
  dayobs<-subset(obs, year==date[[1]][1] & month==date[[1]][2] & day==date[[1]][3])
  
    #get date
    d<-formatC(date[[1]][3],width=2,flag='0')
    m<-formatC(date[[1]][2],width=2,flag='0')
    y<-date[[1]][1]
    
    
    #start and end dates
    obsDate<-as.Date(paste(y,m,d,sep="-"),format="%Y-%m-%d")
    days<-seq(obsDate-14,obsDate+15,by=1)
    lon<-dayobs$lon
    lat<-dayobs$lat
    #function to extract 30 days of data
    climObs<-list()
    for(ii in 1:length(days)){
      climObs[[ii]]<-clim(day=days[ii],lon,lat)
    }
    dat<-as.data.frame(do.call("rbind",climObs))
    dat$obsDate<-obsDate
    dat$ID<-paste(lat,lon)
    climDF<-as.data.frame(summarise(group_by(dat,ID), mean(maxT)))
    climDF$sumPrec<-as.data.frame(summarise(group_by(dat,ID), sum(pre)))[2]$"sum(pre)"
    climDF$meanTmin<-as.data.frame(summarise(group_by(dat,ID), mean(minT)))[2]$"mean(minT)"
    climDF$day<-d
    climDF$month<-m
    climDF$year<-y
    obsClim[[i]]<-climDF
    message(i)
  }

climDF2<-as.data.frame(do.call("rbind",obsClim))

write.csv(climDF2,paste0('/Users/daisy/Google Drive/PhD/BreedingTiming/tables/',
                         'SiloClimateLocationDate30DaySummary', 
                         as.Date(Sys.Date()),'.csv'),row.names=FALSE)
  