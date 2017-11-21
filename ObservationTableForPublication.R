rm(list = ls())

# Breeding quantiles, use to get species
BP<-read.csv("/Users/daisy/Google Drive/PhD/BreedingTiming/tables/BreedingQuantiles2015-10-17.csv")

inc<-read.csv('/Users/daisy/Google Drive/PhD/BreedingTiming/tables/SpeciesOfInterest_2015-10-07.csv')
#keep only species of interest
inc<-subset(inc, remove !=1)
taxon<-inc[c("genus","family","ScientificName")]
BP<-merge(BP,taxon,all.x=TRUE,by.x="Species",by.y="ScientificName")
BP<-subset(BP, !is.na(BP$family))
#unique number of species
species<-unique(BP$Species)


obs.dir<-'/Users/daisy/Google Drive/PhD/Data/Observaitons/Cleaned/Breeding/'

# BP<-merge(BP,taxon,all.x=TRUE,by.x="Species",by.y="ScientificName")
#peak breeding period data
peak<-read.csv(paste0(obs.dir,"PeakBreedingPointOfLayDayOfYear2015-10-17.csv"))
all<-read.csv(paste0(obs.dir,"PointOfLayDayOfYear2015-10-07.csv"))

#make sure only have species interested in

all<-all[all$Scientific.Name %in% species,] # species to keep
all<-all[,c("Scientific.Name","lat","lon","DOY_PL","year","sourceName")]
all$all<-1


#figure out which observations are peak

peak<-peak[peak$Scientific.Name %in% species,] # species to keep
peak<-peak[,c("Scientific.Name","lat","lon","DOY_PL","year","sourceName")]
peak$peak<-1
dat<-merge(all,peak,by = intersect(names(all), names(peak)),all.x=TRUE)


inc2<-subset(inc,remove==1)
RMspecies<-inc2$ScientificName
peak<-peak[peak$Scientific.Name %nin% RMspecies,] #remove unwanted species
BP<-BP[BP$Species %nin% RMspecies,]

# peakSubset <- read.csv("//ad.uws.edu.au/dfshare/HomesHWK$/30022860/Desktop/Daisy/CB/PeakBreedingPointOfLayDayOfYear2015-10-17.csv")
# a<-unique(paste(peakSubset$Scientific.Name,peakSubset$koeppen))
# 
# alldat <- read.csv("/PointOfLayDayOfYear2015-10-07.csv")
# peak$namepeak<-(paste(peak$Scientific.Name,peak$koeppen))
