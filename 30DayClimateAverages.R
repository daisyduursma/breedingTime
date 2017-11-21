rm(list = ls())
library(raster)
library(data.table)
library(dplyr)

silo.dir<-"/Users/daisy/Google Drive/PhD/Data/Spatial/Climate/SILO/"

clim<-c("Precip","Tmax","Tmin")

#Rainfall
for(yr in 1957:2013){
  for (mth in 1:12){
    m<-formatC(mth,width=2,flag='0')
    files<-list.files(paste0(silo.dir,"Daily_Precip/"),
                      pattern=paste0("rai",yr,m),
                      full.names=TRUE)
    if(length(files) <28 | length(files) >31) stop ("incorrect number of daily files")
    message(paste0("year = ",yr," and month = ",mth)) 
    #first file
    r1<-raster(files[1])
    #find total monthly precip
    for (rs in 2:length(files)){
      r1<-r1+raster(files[rs])
    }
    writeRaster(r1,paste0(silo.dir,"Monthly_Precip/rai",yr,m,".tif"))
  }
}

  
#Tmin
for(yr in 1957:2013){
  for (mth in 1:12){
    m<-formatC(mth,width=2,flag='0')
    files<-list.files(paste0(silo.dir,"Daily_Tmin/"),
                      pattern=paste0("Tmin",yr,m),
                      full.names=TRUE)
    if(length(files) <28 | length(files) >31) stop ("incorrect number of daily files")
    message(paste0("year = ",yr," and month = ",mth)) 
    #first file
    r1<-stack(files[1])
    #find total monthly precip
    for (rs in 2:length(files)){
      r1<-stack(r1,raster(files[rs]))
    }
    r1<-mean(r1)
    writeRaster(r1,paste0(silo.dir,"Monthly_Tmin/MMTmin",yr,m,".tif"))
  }
}

#Tmax
for(yr in 1957:2013){
  for (mth in 1:12){
    m<-formatC(mth,width=2,flag='0')
    files<-list.files(paste0(silo.dir,"Daily_Tmax/"),
                      pattern=paste0("Tmax",yr,m),
                      full.names=TRUE)
    if(length(files) <28 | length(files) >31) stop ("incorrect number of daily files")
    message(paste0("year = ",yr," and month = ",mth)) 
    #first file
    r1<-stack(files[1])
    #find total monthly precip
    for (rs in 2:length(files)){
      r1<-stack(r1,raster(files[rs]))
    }
    r1<-mean(r1)
    writeRaster(r1,paste0(silo.dir,"Monthly_Tmax/MMTmax",yr,m,".tif"))
  }
}