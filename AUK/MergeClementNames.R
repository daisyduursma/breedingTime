# 
# traits<-read.csv('/Users/daisy/Google Drive/PhD/BreedingTiming/tables/SpeciesOfInterest_2015-10-07.csv')
# traits<-subset(traits,remove==0)
# clem<-read.csv("~/Google Drive/PhD/Data/Taxonomy/eBird-Clements-2016-integrated-11-Aug.csv", comment.char="#")
# clem<-subset(clem,Category=="species")[,c("sort.v2016","Category","English.name","Scientific.name")]
# 
# traits2<-merge(traits,clem,by.x="Common.me",by.y="English.name",all.x=TRUE)
# traits2$English.name<-traits2$Common.me
# traitsNA<-subset(traits2,is.na(traits2$Scientific.name))[1:43]
# traits2<-subset(traits2,!is.na(traits2$Scientific.name))
# 
# traitsNA<-merge(traitsNA,clem,by.x="ScientificName",by.y="Scientific.name",all.x=TRUE)
# traitsNA$Scientific.name<-traitsNA$ScientificName
# 
# endTraits<-smartbind(traits2,traitsNA)
# 
# write.csv(endTraits,'/Users/daisy/Google Drive/PhD/BreedingTiming/tables/SpeciesOfInterest_2016-10-07.csv')

traits<-read.csv("/Users/daisy/Google Drive/PhD/BreedingTiming/tables/SpeciesOfInterest_2016-10-07.csv")
clem<-read.csv("~/Google Drive/PhD/Data/Taxonomy/eBird-Clements-2016-integrated-11-Aug.csv", comment.char="#")[,c("Scientific.name",        
                                                                                                                  "Order",
                                                                                                                  "Family")]

traits2<-merge(traits,clem,by.x="Clements.Scientific.name",by.y="Scientific.name",all.x=TRUE)
write.csv(traits2,'/Users/daisy/Google Drive/PhD/BreedingTiming/tables/SpeciesOfInterest_2016-10-07.csv',row.names=FALSE)
