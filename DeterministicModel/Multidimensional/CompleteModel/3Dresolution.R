library(deSolve)
library(ggplot2)
library(plyr)

############################################
#### Inj times ########

DACther=F # If we do not consider the DAC therapy

InjEBVTime=c(7,67,127,187,247)*24

if(DACther)
{
  InjDACTime= seq(67,337,by=30)*24
  
}else InjDACTime= NULL

numberDAC <<- 300

numberEBV <<- c(10000,10000)

EBVplaces.injection<-c("px1_py1_pz1","px3_py3_pz3")
DACplaces.injection<-c("px2_py2_pz2","px1_py1_pz3")


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
#p<-c(.4,.2,0.09,0.5,0.1, 3 ,0.15,0.15,0.1)
## Parameters for sick person
p<-c(.4,.2,0.09,0.5,0.15, 1 ,0.1 ,0.1,0.1)

names(p)<- c( "TeE","TrE","Tr2","Te2","TekODC","TrkTe","TekEBV","rec","NKkT")

source("ModelRRMS3D.R")

# yini["EBV_px1_py1_pz1"]=10000
# yini["EBV_px3_py3_pz3"]=10000
# resNoSimm <-lsoda(yini,0:11, funODE, parms=p)

start_time_system = Sys.time()

resNoSimm <- lsoda(yini,seq(from = 0, to = InjTime[1,1] , by = step), funODE, parms=p)

ebv.ind<-1

for(i in 1:length(InjTime[,1]))
{
  yini<-tail(resNoSimm,1)[,-1]
  
  if(InjTime[i,2]=="EBV"){
    yini[paste("EBV_",EBVplaces.injection,sep="")]<-yini[paste("EBV_",EBVplaces.injection,sep="")]+numberEBV
    
  }else  if(InjTime[i,2]=="DAC"){
    yini[dacPos]<-yini[DACplaces.injection]+numberDAC
  }else{
    yini[dacPos]<-yini[DACplaces.injection]+numberDAC
    
    yini[paste("EBV_",EBVplaces.injection,sep="")]<-yini[paste("EBV_",EBVplaces.injection,sep="")]+numberEBV
    
  }
  
  start_time<-InjTime[i,1]
  
  if(i==length(InjTime[,1]))
  {
    end_time<-FinalTime
  }else{
    end_time<-InjTime[(i+1),1]
  }
  #resNoSimm[length(resNoSimm[,2]),1]<- tail(resNoSimm[,1],1)-.2
  resNoSimm <- resNoSimm[-length(resNoSimm[,1]),]
  
  #res <-lsoda(unlist(yini),168:180,funODE,parms=p)
  
  res <-lsoda(yini,seq(from = start_time, to = end_time , by = step),funODE,parms=p)
  
  resNoSimm<-rbind(resNoSimm,res)
}

resNoSimm<-as.data.frame(resNoSimm)
resNoSimm[,"time"]<-resNoSimm[,"time"]/24


end_time = Sys.time()
end_time - start_time_system


