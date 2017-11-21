rm(list = ls())

library(car)
library(doBy)
library(ENmisc)
library(multcomp)
library(circular)
library(MuMIn)
library(lme4)

#working directory
dat.dir<-'/Users/daisy/Google Drive/PhD/Data/Observaitons/Cleaned/Breeding/'
#read in observations
obs<-read.csv(paste0(dat.dir,'PointOfLayDayOfYear2016-09-20.csv'))
#get detials of order
inc<-read.csv('/Users/daisy/Google Drive/PhD/BreedingTiming/tables/SpeciesOfInterest_2016-10-07.csv')
taxon<-inc[c("ScientificName","Order")]

regions<-c(32,13,22,3)

BPregion<-list()
for(rg in 1:length(regions)){#loop through differnt koeppen zones
# #keep only temperate region
 dat<-subset(obs,koeppen==regions[rg])
#get species with more than 100 Egg observations

datEgg<-subset(dat, type == "solitary-egg")
datEgg_sp<-as.data.frame(table(datEgg$Scientific.Name))
datEgg_sp<-as.vector(droplevels(subset(datEgg_sp,Freq>=100)$Var1))
datEgg_sp<-data.frame(Species =datEgg_sp, group ="solitary-egg")
#datEgg<-datEgg[datEgg$Scientific.Name %in% datEgg_sp,]

#get species with more than 100 Unknown observations
#148 species
datUn<-subset(dat, type == "unknown")
datUn_sp<-as.data.frame(table(datUn$Scientific.Name))
datUn_sp<-as.vector(droplevels(subset(datUn_sp,Freq>=100)$Var1))
datUn_sp<-data.frame(Species =datUn_sp, group ="unknown")
datUn<-datUn[datUn$Scientific.Name %in% datUn_sp,]

#get species with more than 100 young observations
#70 species
datYG<-subset(dat, type == "young")
datYG_sp<-as.data.frame(table(datYG$Scientific.Name))
datYG_sp<-as.vector(droplevels(subset(datYG_sp,Freq>=100)$Var1))
datYG_sp<-data.frame(Species =datYG_sp, group ="young")
datYG<-datYG[datYG$Scientific.Name %in% datYG_sp,]

#get species with more than 100 Multivisit observations
#45 species
datmulti<-subset(dat, type == "multi")
datmulti_sp<-as.data.frame(table(datmulti$Scientific.Name))
datmulti_sp<-as.vector(droplevels(subset(datmulti_sp,Freq>=100)$Var1))
if (length(datmulti_sp) !=0){
datmulti_sp<-data.frame(Species =datmulti_sp, group ="multi")
datmulti<-datmulti[datmulti$Scientific.Name %in% datmulti_sp,]

species<-rbind(datYG_sp,datUn_sp,datmulti_sp,datEgg_sp)
} else{
  species<-rbind(datYG_sp,datUn_sp,datEgg_sp)
}
species$region<-regions[rg]

#calculate breeding summaries for differnt data types
datSummary<-list()
for (i in 1:nrow(species)){
      #get species data based on datatype
    spdat<-subset(dat,Scientific.Name==as.character(species[i,]$Species))
    spdat<-subset(spdat,type ==as.character(species[i,]$group) )#get species data
    spdat<-subset(spdat,koeppen ==species[i,]$region )
  
  #calculate summary statistics  
  NoObs<-as.numeric(nrow(spdat))#get observation count
  obsRadians<-spdat$Radians
    #date of quantiles
    quant<-quantile.circular(obsRadians,c(.05,.5,.95),type=8)
    Per5<-round((quant[[1]]*180/pi)/360*365)
    Per50<-round((quant[[2]]*180/pi)/360*365)
    Per95<-round((quant[[3]]*180/pi)/360*365)
    #Breeding Period Length
    StartDate<-as.numeric(Per5)
    EndDate<-as.numeric(Per95)
    BPL<-ifelse (EndDate < StartDate, 365-StartDate+EndDate,EndDate - StartDate )
     # #put the data together and name it
    Perdat<-cbind(species[i,],NoObs,
                  Per5,Per50,Per95,BPL)
    colnames(Perdat)<-c('Species','group','region','Observations','Quantile5',
                        'Quantile50','Quantile95','BreedingPeriod')
    datSummary[[i]]<-Perdat
}

BPregion[[rg]]<-as.data.frame(do.call('rbind', datSummary))

}

datSummary<-do.call("rbind",BPregion)

#add the order
datSummary<-merge(datSummary,taxon,all.x=TRUE,by.x="Species",by.y="ScientificName")


# #make sure values are numeric
# datSummary$BreedingPeriod <-as.numeric(levels(datSummary$BreedingPeriod))[datSummary$BreedingPeriod]
# datSummary$Quantile5 <-as.numeric(levels(datSummary$Quantile5))[datSummary$Quantile5]
# datSummary$Quantile50 <-as.numeric(levels(datSummary$Quantile50))[datSummary$Quantile50]

datSummary$Species<-with(datSummary,paste(Species,region,sep=" "))

#datSummary$Quantile95 <-as.numeric(levels(datSummary$Quantile95))[datSummary$Quantile95]
#order alphebetically
#datSummary <- datSummary[order(datSummary$Species),]

#change peak date so that dates at the beginnig of the year extend over 365
datSummary$Quantile50<-with(datSummary,ifelse(Quantile50<Quantile5,Quantile50+365,Quantile50))

######
modstart<-lmer(Quantile5 ~ region + group + (1|Order),data=datSummary)#add order
Anova(modstart, test.statistic="F")
#find differencs between regions
TukeyRegion1<- glht(modstart, linfct=mcp(group="Tukey"))
summary(TukeyRegion1)
#****unknown significanlty differnt than other data-types

######
modpeak<-lmer(Quantile50 ~ region + group + (1|Order),data=datSummary)#add order
Anova(modpeak, test.statistic="F")
#find differencs between regions
TukeyRegion2<- glht(modpeak, linfct=mcp(group="Tukey"))
summary(TukeyRegion2)
#****unknown significanlty differnt than other data-types

######
modlength<-lmer(BreedingPeriod ~ region + group + (1|Order),data=datSummary)#add order
Anova(modlength, test.statistic="F")
#find differencs between regions
TukeyRegion3<- glht(modlength, linfct=mcp(group="Tukey"))
summary(TukeyRegion3)
#****unknown significanlty differnt than other data-types

summaryBy(BreedingPeriod ~ group, data=datSummary,
          FUN=c(mean,sd)) #find mean and sd




#######################################
#make figure
#######################################

#plot 1 against another with 1 to 1 lines
#fit lm to each ones, add to plot
abline_range <- function(a=NULL,b=NULL,reg=NULL,from=NULL,to=NULL,...){
  
  # Borrowed from abline
  if (!is.null(reg)) a <- reg
  
  if (!is.null(a) && is.list(a)) {
    temp <- as.vector(coefficients(a))
    from <- min(a$model[,2], na.rm=TRUE)
    to <- max(a$model[,2], na.rm=TRUE)
    
    if (length(temp) == 1) {
      a <- 0
      b <- temp
    }
    else {
      a <- temp[1]
      b <- temp[2]
    }
  }
  
  segments(x0=from,x1=to,
           y0=a+from*b,y1=a+to*b,...)
}  


plot.methods <- function(x, y,lmmMod) {
  abline_range(1,1, from=min(subset(x,!is.na(x)&!is.na(y))), to=max(subset(x,!is.na(x)&!is.na(y))), col="black",lty=2)
  coef = fixef(lmmMod)
  abline_range(as.numeric(coef[1]),as.numeric(coef[2]),
               from=min(subset(x,!is.na(x)&!is.na(y))), 
               to=max(subset(x,!is.na(x)&!is.na(y))),col="blue")
  require(MuMIn)
  r.squaredGLMM(lmmMod)[1]#conditional R2, which describes the proportion of variance explained by both the fixed and random factors:
  
  
  legend("bottomright",
         bty="n",
         inset=.1,
         legend=bquote(R^2 == .(round(r.squaredGLMM(lmmMod)[1],2))))
                      
}


plotlabel <- function(txt, where, inset=0.08, font=2, inset.x=inset, inset.y=inset,...){
  u <- par()$usr
  if(grepl("left",where))x <- u[1] + inset.x*(u[2]-u[1])
  if(grepl("right",where))x <- u[2] - inset.x*(u[2]-u[1])
  if(grepl("bottom",where))y <- u[3] + inset.y*(u[4]-u[3])
  if(grepl("top",where))y <- u[4] - inset.y*(u[4]-u[3])
  
  text(x,y,txt,font=font,...)
}



#Start

pdf(file = paste0("/Users/daisy/Google Drive/PhD/BreedingTiming/figures/DataSourcesStartPeakLength",as.Date(Sys.Date()),".pdf"),
    height = 5.5,
    width = 7)



par(mfrow=c(2,3),mar = c(5, 0, 0, 0),oma = c(1, 5, 5, 1))

#peak
peak<-datSummary[,c("Species","Quantile50","group","Order" )]
peak<-reshape(peak, direction="wide", timevar="group", idvar=c("Species","Order"))


plot(peak$Quantile50.young,peak$Quantile50.multi,
     ylab="Multi-visit peak (Julian day)", 
     xlab="Young peak (Julian day)",
     axes=FALSE,frame=TRUE,ylim=c(125,375),xlim=c(100,380))
fit<-lmer(Quantile50.multi~Quantile50.young +(1|Order),data=peak)
plot.methods(peak$Quantile50.young,peak$Quantile50.multi,fit)
Axis(side=2, labels=TRUE,at=seq(150,350,50))
Axis(side=1, labels=TRUE,at=seq(150,350,50))
plotlabel(expression(paste("(", bold("A"), ")")), "topleft")
mtext("Multi-visit peak (Julian day)", side=2, line=3, cex=.7)



plot(peak$"Quantile50.solitary-egg",peak$"Quantile50.multi",
     xlab="Solitary-egg peak (Julian day)",
     axes=FALSE,frame=TRUE,ylim=c(125,375),xlim=c(100,380))
fit2<-lmer(Quantile50.multi~`Quantile50.solitary-egg` +(1|Order),data=peak)
plot.methods(peak$"Quantile50.solitary-egg",peak$Quantile50.multi,fit2)
Axis(side=1, labels=TRUE,at=seq(150,350,50))

plot(peak$"Quantile50.unknown",peak$"Quantile50.multi",
     xlab="Undefined peak (Julian day)",
     axes=FALSE,frame=TRUE,ylim=c(125,375),xlim=c(100,380))
fit3<-lmer(Quantile50.multi~Quantile50.unknown +(1|Order),data=peak)
Anova(fit3)
r.squaredGLMM(fit3)
plot.methods(peak$Quantile50.unknown,peak$Quantile50.multi,fit3)
Axis(side=1, labels=TRUE,at=seq(150,350,50))



#length
lgth<-datSummary[,c("Species","BreedingPeriod","group" )]
lgth<-reshape(lgth, direction="wide", timevar="group", idvar="Species")

plot(lgth$BreedingPeriod.young,lgth$BreedingPeriod.multi,
     ylab="Multi-visit ELP (days)", 
     xlab="Young ELP (days)",
     axes=FALSE,frame=TRUE,ylim=c(75,350),xlim=c(75,350))
fit4<-lmer(lgth$BreedingPeriod.multi~lgth$BreedingPeriod.young +(1|Order),data=peak)
plot.methods(lgth$BreedingPeriod.young,lgth$BreedingPeriod.multi,fit4)
Axis(side=2, labels=TRUE,at=seq(100,300,50))
Axis(side=1, labels=TRUE,at=seq(100,300,50))
plotlabel(expression(paste("(", bold("B"), ")")), "topleft")
mtext("Multi-visit ELP (days)", side=2, line=3, cex=.7)


plot(lgth$"BreedingPeriod.solitary-egg",lgth$"BreedingPeriod.multi",
     xlab="Solitary-egg ELP (days)",
     axes=FALSE,frame=TRUE,ylim=c(75,350),xlim=c(75,350))
fit5<-lmer(lgth$BreedingPeriod.multi~lgth$`BreedingPeriod.solitary-egg` +(1|Order),data=peak)
plot.methods(lgth$"BreedingPeriod.solitary-egg",lgth$"BreedingPeriod.multi",fit5)
Axis(side=1, labels=TRUE,at=seq(100,300,50))

plot(lgth$"BreedingPeriod.unknown",lgth$"BreedingPeriod.multi",
     xlab="Undefined ELP (days)",
     axes=FALSE,frame=TRUE,ylim=c(75,350),xlim=c(75,350))
fit6<-lmer(lgth$BreedingPeriod.multi~lgth$`BreedingPeriod.unknown` +(1|Order),data=peak)
plot.methods(lgth$"BreedingPeriod.unknown",lgth$"BreedingPeriod.multi",fit6)
Axis(side=1, labels=TRUE,at=seq(100,300,50))


dev.off()



