# calculate mean fPAR for each month
library(raster)
dat.dir<-"/Users/daisy/Google Drive/PhD/Data/Spatial/fPAR/Aust_ftot_8km_monthly_1980-2011_v5/"

#make header files - data came without a header file
# 
# hdrfile<-paste0(dat.dir,"Aust_ftot_8km_mth006_198107_v5.hdr")
# 
# files<-list.files(dat.dir, pattern = (".flt"))
# for(i in 1:length(files)){
#   name<-strsplit(files[i],".flt")[[1]][1]
#   file.copy(hdrfile,paste0(dat.dir,name,".hdr"))
# }

#calculate average for each mothb
month<- c("01","02","03","04","05","06","07","08","09","10","11","12")
for(i in 1:length(month)){
  files<-list.files(dat.dir, pattern = paste0(month[i],"_v5.flt"))
  r1<-raster(paste0(dat.dir,files[1]))
  for(ii in 2:length(files)){
    r1<-stack(r1,raster(paste0(dat.dir,files[ii])))
  }
 
  monthly<-mean(r1,na.rm=TRUE)
  
  writeRaster(monthly,paste0("/Users/daisy/Google Drive/PhD/Data/Spatial/fPAR/Aust_ftot_8km_monthly_1980-2011_v5_AVG_month",month[i],".asc"),overwrite=TRUE)
  plot(monthly,main=paste0("avg",month[i]))
}

