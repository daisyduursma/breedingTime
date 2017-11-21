
rm(list = ls())
library(data.table)
# # tables
# inc<-fread('/Users/daisy/Google Drive/PhD/BreedingTiming/tables/SpeciesOfInterest_2015-10-07.csv',data.table=FALSE)
# jetz<-fread('/Users/daisy/Google Drive/PhD/phenology/BirdTreeBLIOCPhyloMasterTax.csv',data.table=FALSE)
# 
# inc2<-merge(inc,jetz,by.x="ScientificName",by.y="Scientific",all.x=TRUE)
# #get complete data
# suninc<-subset(inc2,is.na(inc2$OscSubOsc)==TRUE)
# #get data that does not match
# inc2<-subset(inc2,!is.na(inc2$OscSubOsc))
# #get the data to keep
# fin<-as.data.frame(inc2$"ScientificName") 
# fin$Jetz<-inc2$"ScientificName"
# colnames(fin)<-c("Birdlife","Jetz")
# #merge by common name
# suninc2<-merge(suninc[c(1:43)],jetz,by.x="Common.me",by.y="English",all.x=TRUE)
# fin2<-suninc2[c("ScientificName","Scientific")]
# colnames(fin2)<-c("Birdlife","Jetz")
# 
# enddat<-rbind(fin,fin2)
# write.csv(enddat,'/Users/daisy/Google Drive/PhD/phenology/JetzBirdlifeNames.csv',row.names=FALSE)



########################
#read in data created above
########################
jetz<-fread('/Users/daisy/Google Drive/PhD/phenology/JetzBirdlifeNames.csv',data.table=FALSE)


########################
#get names for temperate region
########################

temperate<-fread("/Users/daisy/Google Drive/PhD/BreedingTiming/tables/BreedingQuantiles2015-10-17.csv",data.table=FALSE)
temperate<-as.data.frame(subset(temperate,Region=="Temperate"))
temperate$allBreeding<-1
temperate<-temperate[c("Species","allBreeding")]

jetz<-merge(jetz,temperate,all.x=TRUE,by.x="Birdlife",by.y="Species")



########################
#get names for temperate region
########################


#get names for el Nino/la Nina analyses
  #working directory
  dat.dir<-'/Users/daisy/Google Drive/PhD/Data/Observaitons/Cleaned/Breeding/'
  #read in observations
  dat<-read.csv(paste0(dat.dir,'PointOfLayDayOfYear2015-10-07.csv'))
  #add koeppen zones
  koeppen<-raster(paste0('/Users/daisy/Google Drive/PhD/Data/Spatial/BOM_climate_zones/kpngrp_major/koepenReclassified.asc'),
                  proj4string=CRS("+proj=longlat +datum=WGS84 +ellps=WGS84"))
  #combine equatorial and tropical: 41 Equatorial, 35 Tropical, extract koeppen zone
  m <- c(35, 41, 35)
  rclmat <- matrix(m, ncol=3, byrow=TRUE)
  koeppen<- reclassify(koeppen, rclmat)
  dat$koeppen<-extract(koeppen,data.frame(cbind(dat$lon, dat$lat)))
  #convert data to radians to make circular vector
  dat$Radians <-(dat$DOY_PL/366*360)*pi / 180
  #traits for species
  inc<-read.csv('/Users/daisy/Google Drive/PhD/BreedingTiming/tables/SpeciesOfInterest_2015-10-07.csv')
  #get list of species not to include
  inc2<-subset(inc,remove==1)
  RMspecies<-inc2$ScientificName
  dat<-dat[dat$Scientific.Name %nin% RMspecies,]
  inc<-inc[inc$ScientificName %nin% RMspecies,]#remove unwanted species
  #Koeppen zones: 35-Tropical, 32-Subtropical, 22-Desert,3-Grassland,3-Temperate
  TEMP<-subset(dat, koeppen ==3) #temperate region Breeding Period
  TEMP$month<-formatC(TEMP$month,width=2,flag='0') #format months 1 becomes 01
  TEMP$date<-paste0(TEMP$year,TEMP$month) #make date string that matches NOAA SOI
  TEMP<-subset(TEMP, date >=195101) #limit time to same as NOAA data
  #clades<-subset(inc,select=c("ScientificName","Patch.clade"))
  #TEMP<-merge(TEMP,clades, by.x ="Scientific.Name", by.y ="ScientificName")
  
  #read in Southern Oscillation Index
  #SOI<-read.csv("/Users/daisy/Google Drive/PhD/Data/SOI/NOAA_SouthernOscillationIndex.csv")
  SOI<-read.csv("/Users/daisy/Google Drive/PhD/Data/SOI/NOAA_OceanicNinoIndex5MonthPhase.csv")
  colnames(SOI)<-c("Year",1:12)
  SOI<-melt(SOI, id="Year")
  SOI$month<-formatC(SOI$variable,width=2,flag='0') #format months 1 becomes 01
  SOI$date<-paste0(SOI$Year,SOI$month) #make date string that matches NOAA SOI
  TEMP<-merge(TEMP,SOI, by.x="date",by.y="date",all.x=TRUE)
  
  
  ##########################################
  ############ step 4. get species with more than 50 obs in each phase ########
  ##########################################
  
  #get el Nino obs and species with more than 200 observations
  eN<-subset(TEMP,value=='E')
  spSumeN<-as.data.frame(table(eN$Scientific.Name))
  sp_eN<-as.vector(droplevels(subset(spSumeN,Freq>=100)$Var1))
  lN<-subset(TEMP,value=='L')
  spSumlN<-as.data.frame(table(lN$Scientific.Name))
  sp_lN<-as.vector(droplevels(subset(spSumlN,Freq>=100)$Var1))
  nN<-subset(TEMP,value=='N')
  spSumnN<-as.data.frame(table(nN$Scientific.Name))
  sp_nN<-as.vector(droplevels(subset(spSumnN,Freq>=100)$Var1))
  #name of species in both, limit data to these species
  fin_sp<-as.data.frame(intersect(sp_eN,sp_lN))
fin_sp$NinaNino<-1
colnames(fin_sp)<-c("Species","NinaNino")

jetz<-merge(jetz,fin_sp,all.x=TRUE,by.x="Birdlife",by.y="Species")


write.csv(jetz,'/Users/daisy/Google Drive/PhD/phenology/JetzBirdlifeNames.csv',row.names=FALSE)
