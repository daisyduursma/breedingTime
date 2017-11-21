


#clean and combine breeding observation data from all museums, ATLAS, Nest Record Scheme
#seperate epoch dat into nest building, egg, young, and unknown event
#fix lat. values that are not negative
#check observations falll within Australia based on eMAST data
#calculate epoch date since 1970-1-1
#write table of species with more than 100 breeding observations


rm(list = ls())

# library(raster)
library(car) 
library(stringr)
library(Hmisc)
library(lubridate)
library(gtools)

indir<-'/Users/daisy/Google Drive/PhD/Data/Observaitons/Raw/egg'


abbbsA<-read.csv(paste0(indir,"/ABBBS/pullusBandingRecordsPartA_dd.csv"))
abbbsB<-read.csv(paste0(indir,"/ABBBS/pullusBandingRecordsPartB_dd.csv"))
abbbs<-rbind(abbbsA,abbbsB)
abbbs<-abbbs[!duplicated(abbbs[,c("SCIENTIFIC_NAME","Day","Month","Year","LAT","LON","BANDER")]),]
abbbs$sourceName<-'ABBBS'
abbbs$epoch<-as.numeric(as.Date(paste0(abbbs$Year,'-',abbbs$Month,'-',abbbs$Day)))
abbbs$type<-with(abbbs, ifelse (AGE == "NESTLING","youngNest",'unknown'))
abbbs$startUnknown<-with (abbbs, ifelse (type=='unknown', abbbs$epoch,NA))  
abbbs$startYoung<-with (abbbs, ifelse (type=='startYoung', abbbs$epoch,NA))
abbbs$lat<-abbbs$LAT
abbbs$lon<-abbbs$LON
abbbs$Scientific.Name<-abbbs$SCIENTIFIC_NAME 
abbbs$sourceName <- 'ABBBS'

abbbs<-abbbs[,c('Scientific.Name','lat','lon', 'sourceName', 'type' ,'startYoung','startUnknown','BANDER')]
eggdat<-abbbs

# returns string w/o leading or trailing whitespace
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

#finish cleaning egdat
eggdat$Scientific.Name<- capitalize(trim(gsub('  ', ' ', eggdat$Scientific.Name,
                                              fixed=TRUE)))#check for extra spaces and capatlize
eggdat<-subset(eggdat, !is.na(lat) & !is.na(lon) & !is.na(Scientific.Name))#make sure lat and longs are given

#make sure all latitudes are negative
latfix<-subset(eggdat, lat > 0)
latfix$lat<-latfix$lat*-1
latgood<-subset(eggdat, lat <= 0)
eggdat<-rbind(latfix, latgood)

#remove duplicates
eggdat<-eggdat[!duplicated(eggdat),]

#make sure years are as expected
syn<-subset(read.csv('/Users/daisy/Google Drive/PhD/Data/Observaitons/Raw/namesResolved.csv'), use ==1)
bad<-subset(read.csv('/Users/daisy/Google Drive/PhD/Data/Observaitons/Raw/namesResolved.csv'), use ==0)
eggdat<-eggdat[eggdat$Scientific.Name %nin% bad$scn1,]
fix<-eggdat[eggdat$Scientific.Name %in% syn$scn1,]
nofix<-eggdat[eggdat$Scientific.Name %nin% syn$scn1,]#[c(1:8)]
fix<-merge(fix,syn,by.x = 'Scientific.Name',by.y='scn1')
fix$Scientific.Name <-fix$garnett_name


fix<-fix[,c('Scientific.Name','lat','lon', 'sourceName', 'type' ,'startYoung','startUnknown','BANDER')]
eggdat<-smartbind(fix,nofix)
eggdat<-eggdat[!duplicated(eggdat),]

#make sure data is within continental Australia
aus<-raster('/Users/daisy/Google Drive/PhD/Data/Spatial/Climate/EMASTBiovars/bio_1.asc',
            crs = '+proj=longlat +datum=WGS84')
xy<-cbind(eggdat$lon,eggdat$lat)
eggdat$outsideAustralia<-is.na(extract(aus,xy))
eggdat<-subset(eggdat,outsideAustralia==FALSE)

#see if they are species of interest

inc<-read.csv('/Users/daisy/Google Drive/PhD/BreedingTiming/tables/SpeciesOfInterest_2015-10-07.csv')
#make sure only have species interested in
inc2<-subset(inc,remove==1)
RMspecies<-inc2$ScientificName
eggdat<-eggdat[eggdat$Scientific.Name %nin% RMspecies,] #remove unwanted species

a<- trim(gsub('  ', ' ', eggdat$BANDER,fixed=TRUE))#remove white spaces
a<-as.data.frame(unique(a))#list of names need to include

write.csv(a,"~/Google Drive/PhD/BreedingTiming/tables/ABBBSBanders.csv")
