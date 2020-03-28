library(deSolve)
library(ggplot2)
library(plyr)

############################################
#### Inj times ########

InjEBVTime=c(7,67,127,187,247)*24
InjDACTime= NULL

### Let define the number of DAC and EBV cells injected each time

numberDAC <<- 0
numberEBV <<- 1000

###############


###### structure used for the stopping and restarting lsoda

L<-length(c(InjDACTime,InjEBVTime))

InjTime<-data.frame("Times" = numeric(L), "kind" =character(L) )
InjTime$Times = c(InjDACTime,InjEBVTime)
InjTime$kind=rep(c("DAC","EBV"),c(length(InjDACTime),length(InjEBVTime)) )
InjTime<-InjTime[with(InjTime, order(Times)),]
rep<-count(InjTime[,1])

if( TRUE%in%c(rep[,2]>1) )
{
  rep_pos<-which(rep[rep[,2]>1,1]==InjTime)
  InjTime$kind[rep_pos]="BOTH"
  InjTime<-unique(InjTime)
}

##################################################

FinalTime<-365*24

step=24

## Parameters for healthy person
#p<-c(.4,.2,0.09,0.5,0.1, 3 ,0.15,0.1,0.1)
## Parameters for sick person
p<-c(.4,.2,0.09,0.5,0.15, 1 ,0.1 ,0.1,0.1)

names(p)<- c( "TeE","TrE","Tr2","Te2","TekODC","TrkTe","TekEBV","rec","NKkT")


################################################################
########## Model reading and resolution

source("ModelRRMS.R")

res1 <-lsoda(yini,seq(from = 0, to = InjTime[1,1], by = step), funODE, parms=p)

for(i in 1:length(InjTime[,1]))
{
  yini<-tail(res1,1)[,-1]
  
  if(InjTime[i,2]=="EBV"){
    ebvPosition<-sample(.GlobalEnv$ebvPos,1)
    yini[ebvPosition]<-yini[ebvPosition]+numberEBV
  }else  if(InjTime[i,2]=="DAC"){
    yini[dacPos]<-yini[.GlobalEnv$dacPos]+numberDAC
  }else{
    yini[dacPos]<-yini[.GlobalEnv$dacPos]+numberDAC
    ebvPosition<-sample(.GlobalEnv$ebvPos,1)
    yini[ebvPosition]<-yini[ebvPosition]+numberEBV
  }
  
  start_time<-InjTime[i,1]
  
  if(i==length(InjTime[,1]))
  {
    end_time<-FinalTime
  }else{
    end_time<-InjTime[(i+1),1]
  }
  #res1[length(res1[,2]),1]<- tail(res1[,1],1)-.2
  res1 <- res1[-length(res1[,1]),]
  res <-lsoda(yini,seq(from = start_time, to = end_time, by = step),funODE,parms=p)
  
  res1<-rbind(res1,res)
}

res1<-as.data.frame(res1)
res1[,"time"]<-res1[,"time"]/24

ggplot(res1,aes(x=time))+geom_line(aes(y=res1[,ebvPos+1]))
ggplot(res1,aes(x=time))+geom_line(aes(y=res1[,nkPos+1]))
ggplot(res1,aes(x=time))+geom_line(aes(y=res1[,"ODC_le1"]))

ggplot(res1,aes(x=time))+geom_line(aes(y=res1[,teffPos+1],col="Teff") )+geom_line(aes(y=res1[,tregPos+1],col="Treg"))

save(res1, "./Desktop/CMISF2019/prova/",)
