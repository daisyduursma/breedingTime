
rm(list = ls())
# library(raster)
# library(circular)
# library(gtools)
# library(Hmisc)
# library(plyr)
# library(Hmisc)
# library(multcomp)
# library(gplots)
# library(reshape2)
# library(reporttools)
# #figure
# library(lme4)
# #library(afex)
library(ENmisc)
library(doBy)
# library(arm)
# library(car)




##########################################
############ step 12. Figure for Species on the move- Breeding Period by Phase   ########
##########################################


boxSlopes <- function(condition,group,labelY,sig,sigsize,figNumber,Sp,dat,axesTF,at,labels){
  wtd.boxplot(dat[,paste(condition)]~dat[,paste(group)],
              ylab=labelY, 
              main="",names=c("","","","",""),outcol="white", 
              boxwex=0.5,border = c("black","black"),
              cex.lab=1, axes=axesTF)
  box()
  axis(1, at=1:5,labels=FALSE)
  mtext("Dst", side=1, line=1, at=1, cex=1)
  mtext("Grs", side=1, line=1, at=2, cex=1)
  mtext("Tmp", side=1, line=1, at=3, cex=1)
  mtext("Sbt", side=1, line=1, at=4, cex=1)
  mtext("Trp", side=1, line=1, at=5, cex=1)
  if(axesTF==FALSE){
    axis(2, at=at, 
         labels=labels)
  }
  mtext(figNumber, side=3, line=1, at=0,cex=1)
  
  
  #add lines for species in Dst and Grs
  
  for(s in 1:length(Sp)){
    D<-subset(dat,Region == "Desert" & Species == Sp[s])[,paste(condition)]
    G<-subset(dat,Region == "Grassland" & Species == Sp[s])[,paste(condition)]
    if (length(D) == 0 | length(G) == 0){
      next
    } else {
      color=ifelse(subset(dat,Region == "Desert" & Species == Sp[s])[,"DoD"]=="Precocial", "dodgerblue3", "grey60")
      segments(1,D,2,G, col=color, lwd=0.8, lty=1)
    }
  }
  
  for(s in 1:length(Sp)){
    G<-subset(dat,Region == "Grassland" & Species == Sp[s])[,paste(condition)]
    Tm<-subset(dat,Region == "Temperate" & Species == Sp[s])[,paste(condition)]
    if (length(G) == 0 | length(Tm) == 0){
      next
    } else {
      color=ifelse(subset(dat,Region == "Grassland" & Species == Sp[s])[,"DoD"]=="Precocial", "dodgerblue3", "grey60")
      segments(2,G,3,Tm, col=color, lwd=0.8, lty=1)
    }
  }
  
  for(s in 1:length(Sp)){
    S<-subset(dat,Region == "Subtropical" & Species == Sp[s])[,paste(condition)]
    Tm<-subset(dat,Region == "Temperate" & Species == Sp[s])[,paste(condition)]
    if (length(S) == 0 | length(Tm) == 0){
      next
    } else {
      color=ifelse(subset(dat,Region == "Subtropical" & Species == Sp[s])[,"DoD"]=="Precocial", "dodgerblue3", "grey60")
      segments(3,Tm,4,S, col=color, lwd=0.8, lty=1)
    }
  }

  for(s in 1:length(Sp)){
    S<-subset(dat,Region == "Subtropical" & Species == Sp[s])[,paste(condition)]
    Tp<-subset(dat,Region == "Tropical" & Species == Sp[s])[,paste(condition)]
    if (length(S) == 0 | length(Tp) == 0){
      next
    } else {
      color=ifelse(subset(dat,Region == "Subtropical" & Species == Sp[s])[,"DoD"]=="Precocial", "dodgerblue3", "grey60")
      segments(4,S,5,Tp, col=color, lwd=0.8, lty=1)
    }
  }
  #add in the box plots again so that they are not covered by the lines
  wtd.boxplot(dat[,paste(condition)]~dat[,paste(group)],
              ylab=labelY, 
              main="",names=c("","","","",""),outcol="white", 
              boxwex=0.5,border = c("black","black"),
              cex.lab=1, axes=axesTF, add = TRUE, col ="#FFFFFF00" )
  mtext(sig, side=3, line=-2.5, at=3,  adj=0.5,cex=sigsize)
}






############################

dat<-read.csv("/Users/daisy/Google Drive/PhD/BreedingTiming/manuscript/TablesFigures/SupplementaryData1.csv")
dat<-subset(dat, No..records>=100)

#prep data for those species that occur in more than 3 biomes
subDat<-list()
for (i in 1:length(sp)){
  spdat<-subset(dat,Species==sp[i])
  if(nrow(spdat) <= 2) {
    next
  } else {
    subDat[[i]]<-spdat
  
  }
}
dat<-as.data.frame(do.call('rbind', subDat))
#make grouping code for plots
des<-subset(dat,Region=="Desert")
des$biome<-1
grs<-subset(dat,Region=="Grassland")
grs$biome<-2
sbt<-subset(dat,Region=="Subtropical")
sbt$biome<-4
Tmp<-subset(dat,Region=="Temperate")
Tmp$biome<-3
Trp<-subset(dat,Region=="Tropical")
Trp$biome<-5
dat<-rbind(des,grs,sbt,Tmp,Trp) 

#add precocial and altrical
inc<-read.csv('/Users/daisy/Google Drive/PhD/BreedingTiming/tables/SpeciesOfInterest_2015-10-07.csv')[,c("Common.me",
                                                                                                      "family",
                                                                                                      "order",
                                                                                                      "DoD")]
dat<-merge(dat,inc,by.x="Common.name",by.y="Common.me",all.x=TRUE)







#######################################


#set up PDF
pdf(file = "/Users/daisy/Google Drive/PhD/BreedingTiming/figures/SpagattieGraphIntraSpecificVariation.pdf",
    width = 4, height = 9)

par(mfrow=c(3,1))


#Three figures 1 row 2 columns, breeding period, start
boxSlopes(condition="Egg.laying.period",
          group="biome",
          labelY="Length (days)",
          sig="*",
          sigsize=1.5,
          figNumber="(a)",
          Sp=unlist(unique(dat$Species)),
          dat=dat,
          axesTF=TRUE)

 boxSlopes(condition="X5th.percentile",
           group="biome",
           labelY="Start of ELP",
           sig="n.s.",
           sigsize=.9,
           figNumber="b",
           Sp=unlist(unique(dat$Species)),
           dat=dat,
           axesTF=FALSE,
           at=c(1,32,61,92,121,152,182,214,245,275,335),
           labels=c("Jan", "Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov"))

 dat$MODIFIED_Quantile95<-with(dat,ifelse(X95th.percentile<200,X95th.percentile+365,X95th.percentile))
 boxSlopes(condition="MODIFIED_Quantile95",
           group="biome",
           labelY="Termination of ELP",
           sig="*",
           sigsize=1.5,
           figNumber="(c)",
           Sp=unlist(unique(dat$Species)),
           dat=dat,
           axesTF=FALSE,
           at=c(213,245,275,305,335,366,397,425,456,487,518,548),
           labels=c("Aug","Sep","Oct","Nov","Dec","Jan","Feb","Mar","Apr","May","Jun","July"))


 dev.off()


##########################
 
 #How many species have shorter breeding periods on the temperate region than any other region
 
 #######################
 
 pre<-subset(dat, DoD == "Precocial")
 SE <- function(x)sd(x)/sqrt(length(x))
 mn <- summaryBy(Egg.laying.period ~ Region, data=pre,
                 FUN=c(mean,SE,sd,min,max))
 minPre<-list()
 for(ii in 1: length(unlist(unique(pre$Species)))){
   subPre<-subset(pre,Species==unlist(unique(pre$Species))[ii])
   minPre[[ii]]<-subset(subPre,Egg.laying.period==min(subPre$Egg.laying.period))
 }
 do.call("rbind",minPre)
 
 
 #altricial
 alt<-subset(dat, DoD == "Altricial")
 SE <- function(x)sd(x)/sqrt(length(x))
 summaryBy(Egg.laying.period ~ Region, data=alt,
                 FUN=c(mean,SE,sd,min,max))
 
 minAlt<-list()
 for(ii in 1: length(unlist(unique(alt$Species)))){
   subAlt<-subset(alt,Species==unlist(unique(alt$Species))[ii])
   short<-subset(subAlt,Egg.laying.period==min(subAlt$Egg.laying.period))
   Tem<-subset(subAlt,Region=="Temperate")[,"Egg.laying.period"]
   short$diff<-short$Egg.laying.period-Tem
   minAlt[[ii]]<-short
 }
 minAlt<-do.call("rbind",minAlt)
 
 table(minAlt$Region)
 
 nrow(subset(minAlt,diff <=-20))
 
 
 ############################################
 
 #where do species have the earliest start of their breeding period
 
 ############################################
 
 pre<-subset(dat, DoD == "Precocial")
 summaryBy(X5th.percentile ~ Region, data=pre,
                 FUN=c(mean,SE,sd,min,max))
 minPre<-list()
 for(ii in 1: length(unlist(unique(pre$Species)))){
   subPre<-subset(pre,Species==unlist(unique(pre$Species))[ii])
   minPre[[ii]]<-subset(subPre,X5th.percentile==max(subPre$X5th.percentile))
 }
 do.call("rbind",minPre)
 
 
 #altricial
 alt<-subset(dat, DoD == "Altricial")
 summaryBy(X5th.percentile ~ Region, data=alt,
           FUN=c(mean,SE,sd,min,max))
 
 
 earlAlt<-list()
 lateAlt<-list
 for(ii in 1: length(unlist(unique(alt$Species)))){
   subAlt<-subset(alt,Species==unlist(unique(alt$Species))[ii])
   earliest<-subset(subAlt,X5th.percentile==min(subAlt$X5th.percentile))
   latest<-subset(subAlt,X5th.percentile==max(subAlt$X5th.percentile))
   earlAlt[[ii]]<-earliest
   lateAlt[[ii]]<-latest
 }
 earlAlt<-do.call("rbind",earlAlt)
 lateAlt<-do.call("rbind",lateAlt)
 
 table(earlAlt$Region)
 table(lateAlt$Region)
 
 nrow(subset(minAlt,diff <= 20 & diff >0))
 




