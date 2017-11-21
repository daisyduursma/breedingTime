# calculate mean GPP for each month
library(raster)
dat.dir<-"/Users/daisy/Google Drive/PhD/Data/Spatial/GPP/Original"

#calculate average for each month
month<- c("01","02","03","04","05","06","07","08","09","10","11","12")
for(i in 5:length(month)){
  files<-list.files(dat.dir, pattern = paste0(month[i],".nc"))
  try(if(length(files)!=11) stop("unexpected number of files"))
  r1<-raster(paste0(dat.dir,"/",files[1]),varname = 'total')
  for(ii in 2:length(files)){
    r1<-stack(r1,raster(paste0(dat.dir,"/",files[ii]),varname = 'total'))
  }
 
  monthly<-mean(r1,na.rm=TRUE)
  writeRaster(monthly,paste0("/Users/daisy/Google Drive/PhD/Data/Spatial/GPP/AverageMonthly/Hume_GPP_250m_v5_monthlyAverage_",month[i],".nc"),overwrite=TRUE)
}

