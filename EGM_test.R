January, 0; February, 1; March, 285; April, 1637; May, 1413; June, 726;
July, 163; August, 56; September, 4; October, 0; November, 0; December, 0;
This content downloaded from 137.111.13.200 on Fri, 31 Jul 2015 00:23:33 UTC
All use subject to JSTOR Terms and Conditions
BIRD SPECIES DIVERSITY 395
for a total of 4285. The month diversity is 1.434 and el.434 = 4.19, so that
England has the equivalent of 4.19 equally good months for nesting.
In order to find how many shifts there are in the season for nesting

=4.19
# 
# CalcEGM <- function(mothlyObs){
#   EGM<-exp( -( sum((mothlyObs/sum(mothlyObs))*log(mothlyObs/sum(mothlyObs))) ) )
#  return(EGM)
#}

CalcEGM <- function(monthlyObs){
  a<-monthlyObs/sum(monthlyObs)
  b<- log(a)
  b[is.infinite(b)] <- 0 
  EGM<-exp( -( sum(a*b) ) )
  return(EGM)
}
dat1<-c(0,1,285,1637,1413,726,163,56,4,0,0,0)

#expected EGM =6.78
dat2<-c(28,68,155,401,344,148,91,68,15,11,7,21)
#mothlyObs<-c(1,285,1637,1413,726,163,56)

CalcEGM(dat1)
CalcEGM(dat2)

