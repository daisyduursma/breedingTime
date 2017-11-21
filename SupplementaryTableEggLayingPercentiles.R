#get data ready for supplementary information, make sure names are up to date. 

rm(list = ls())

# Breeding quantiles
BP<-read.csv("/Users/daisy/Google Drive/PhD/BreedingTiming/tables/BreedingQuantiles2015-10-17.csv")

inc<-read.csv('/Users/daisy/Google Drive/PhD/BreedingTiming/tables/SpeciesOfInterest_2015-10-07.csv')
#keep only species of interest
inc<-subset(inc, remove !=1)
taxon<-inc[c("genus","family","ScientificName")]
BP<-merge(BP,taxon,all.x=TRUE,by.x="Species",by.y="ScientificName")
BP<-subset(BP, !is.na(BP$family))
#unique number of species
length(unique(BP$Species))
#check names are okay with Australian Bird Data Version 1
birdData<-read.csv("/Users/daisy/Google Drive/PhD/birdTraits/PublishedData/Australian_Bird_Data_Version_1.csv")
#keep only full species
birdData2<-subset(birdData, is.na(X6_Subspecies_name_2))[c(1:5)]

birdData2$Species<-with(birdData2,paste(X4_Genus_name_2,X5_Species_name_2,sep=" "))
BP2<-merge(BP,birdData2,all.x=TRUE,by="Species")
bad<-subset(BP2,is.na(X4_Genus_name_2))
good<-subset(BP2,!is.na(X4_Genus_name_2))
good<-good[,c("Species",
            "CommonName",
            "Region",
            "Quantile5",
            "Quantile50",
            "Quantile95",
            "BreedingPeriod",
            "ObservationCount")]

i <- sapply(good, is.factor)
good[i] <- lapply(good[i], as.character)


bad<-merge(bad,birdData2,all.x=TRUE,by.x="CommonName",by.y="X3_Taxon_common_name_2")
bad$Species<-bad$Species.y
bad<-bad[,c("Species",
            "CommonName",
            "Region",
            "Quantile5",
            "Quantile50",
            "Quantile95",
            "BreedingPeriod",
            "ObservationCount")]
#get rid of factors
i <- sapply(bad, is.factor)
bad[i] <- lapply(bad[i], as.character)

all<-rbind(good,bad)

write.csv(all,"/Users/daisy/Google Drive/PhD/BreedingTiming/manuscript/TablesFigures/SupplementaryData1.csv")

