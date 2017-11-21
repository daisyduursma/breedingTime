
EGM <- read.csv("//ad.uws.edu.au/dfshare/HomesHWK$/30022860/Desktop/Daisy/CB/EGM_breeding2016-03-04.csv")

#   35 Tropical and Equatorial
#   32 Subtropical
#   22 Desert
#   13 Grassland
#   3 Temperate

s1 <- read.csv("//ad.uws.edu.au/dfshare/HomesHWK$/30022860/Desktop/Daisy/CB/Table S1.csv")

df<-merge(s1,EGM,all.x=TRUE,by.x=c("Species","Region"), by.y=c("Species","Biome"))
bad<-subset(df,is.na(Family.scientific.name))
good<-subset(df,!is.na(Family.scientific.name))

i <- sapply(good, is.factor)
good[i] <- lapply(good[i], as.character)

good<-good[,c("Species",
            "Common.name",
            "Region",
            "X5th.percentile",         
            "X50th.percentile",         
            "X95th.percentile", 
            "Egg.laying.period",
            "No..records",
            "EGM")]

colnames(good)<-c("Species",
                  "Common.name",
                  "Region",
                  "X5th.percentile",         
                  "X50th.percentile",         
                  "X95th.percentile", 
                  "Egg.laying.period",
                  "No.records",
                  "EGM")


bad<-merge(bad,EGM,all.x=TRUE,by.x=c("Common.name","Region"),by.y=c("Taxon.common.name","Biome"))
bad<-bad[,c("Species.x",
            "Common.name",
            "Region",
            "X5th.percentile",         
            "X50th.percentile",         
            "X95th.percentile", 
            "Egg.laying.period",
            "No..records",
            "EGM.y")]
#get rid of factors
i <- sapply(bad, is.factor)
bad[i] <- lapply(bad[i], as.character)

colnames(bad)<-c("Species",
                  "Common.name",
                  "Region",
                  "X5th.percentile",         
                  "X50th.percentile",         
                  "X95th.percentile", 
                  "Egg.laying.period",
                  "No.records",
                  "EGM")


all<-rbind(good,bad)

write.csv(all,"//ad.uws.edu.au/dfshare/HomesHWK$/30022860/Desktop/Daisy/CB/TableS120160307.csv")

