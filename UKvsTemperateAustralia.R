rm(list = ls())
library(Hmisc)
library(doBy)

# Simple function for placing labels on a figure.
plotlabel <- function(txt, where, inset=0.08, font=2, inset.x=inset, inset.y=inset,...){
  u <- par()$usr
  if(grepl("left",where))x <- u[1] + inset.x*(u[2]-u[1])
  if(grepl("right",where))x <- u[2] - inset.x*(u[2]-u[1])
  if(grepl("bottom",where))y <- u[3] + inset.y*(u[4]-u[3])
  if(grepl("top",where))y <- u[4] - inset.y*(u[4]-u[3])
  
  text(x,y,txt,font=font,...)
}


#Australia data
BP<-read.csv("/Users/daisy/Google Drive/PhD/BreedingTiming/tables/BreedingQuantiles2015-10-17.csv")
inc<-read.csv('/Users/daisy/Google Drive/PhD/BreedingTiming/tables/SpeciesOfInterest_2015-10-07.csv')
#get list of species not to include
inc2<-subset(inc,remove==1)
RMspecies<-inc2$ScientificName
BP<-BP[BP$Species %nin% RMspecies,]

#################
# What are are the differences in the NHTR vs. Temperate Australia? 
#################

TEMPdat<- subset(BP, Region=="Temperate",select=c("BreedingPeriod","Region")) #AUS temp. breeding period length
UK<-subset(read.csv('/Users/daisy/Google Drive/PhD/BreedingTiming/tables/Joys_etal_Breeding_periods_England.csv')
             ,select="length.of.breeding") #UK breeding period length
colnames(UK)<-"BreedingPeriod"
UK$Region<-"UK"
df<-droplevels(rbind(TEMPdat,UK))#dataframe of UK and AUS temp

#Run ANOVA
bp<-lm(formula = BreedingPeriod ~ Region, data = df,REML=T)
summary(bp)
Anova(bp)

#means and SD of regions
summaryBy(BreedingPeriod ~ Region, data=df,
                  FUN=c(mean,sd))
 


#################
# #transparent histograms of breeding period length, 1 column wide
#################
pdf(file = "/Users/daisy/Google Drive/PhD/BreedingTiming/figures/PeakBreedinPeriodUKandTemperate201604128.pdf",
    height = 4,
    width = 3.5)

par(mfrow=c(2,1), oma=c(4,4,0,0), mar=c(1,1,1,1), las=1, cex=.7, xaxs="i",
    yaxs = "i")
  #set the intervals for the histograms
  breaks <- seq(0,350,25)
  #PREP DATA
  d <- split(df, df$Region)
  a = d[[1]][,1]
  b = d[[2]][,1]
  #UK
  hist(b,xlim=c(0,350), 
       col=rgb(0, 0, 0,0.5),
       breaks = breaks,
       main=NULL, 
       axes=FALSE,
       ylim=c(0,70))
  plotlabel("(a)", "topright")
  Axis(side=1, labels=FALSE,at=seq(0,350,25))
  Axis(side=2, at=c(0,10,20,30,40,50,60,70),labels=c(0,10,20,30,40,50,60,70))
  #AUS temperate
  hist(a,xlim=c(0,350), 
       col=rgb(0, 0, 0,0.5),
       breaks = breaks,
       main=NULL,
       axes=FALSE,
       ylim=c(0,70))
  plotlabel("(b)", "topright")
  Axis(side=1, labels=FALSE,at=seq(0,350,25))
  Axis(side=1, labels=seq(0,350,50),at=seq(0,350,50))
  Axis(side=2, at=c(0,10,20,30,40,50,60,70),labels=c(0,10,20,30,40,50,60,70))
  #y and y labels
  mtext("Frequency", 2, 2, outer=TRUE,las=0)
  mtext("Breeding period (days)", 1, 2, outer=TRUE)
  
  dev.off()
  
  
  
  
  
  