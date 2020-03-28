library(deSolve)
library(ggplot2)
library(plyr)

### If DACther = True then model with dac therapy

DACther=F

############################################
#### Inj times ########

InjEBVTime=c(7,67,127,187,247)*24
if(DACther)
{
  InjDACTime= seq(60,360,by=30)*24
  
  }else InjDACTime= NULL

### Let define the number of DAC and EBV cells injected each time

numberDAC <<- 30
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

y_names<-c("EBV","Teff","Treg","ODC_le1","ODC_le2","ODC_le3","ODC_le4","ODC_le5","NK","IL2","DAC","Resting_Teff","Resting_Treg","EffectorMemory","Resting_Treg_temp","Resting_Teff_temp","NK_temp")

names(yini)= y_names


dacPos<-c( which(y_names %in% grep("DAC", y_names, value=T)) )
ebvPosition<-c( which(y_names %in% grep("EBV", y_names, value=T)) )
###########
res1 <-lsoda(yini,seq(from = 0, to = InjTime[1,1], by = step), funODE, parms=p)

for(i in 1:length(InjTime[,1]))
{
  yini<-tail(res1,1)[,-1]
  
  if(InjTime[i,2]=="EBV"){
   
    yini[ebvPosition]<-yini[ebvPosition]+numberEBV
  }else  if(InjTime[i,2]=="DAC"){
    yini[dacPos]<-yini[.GlobalEnv$dacPos]+numberDAC
  }else{
    yini[dacPos]<-yini[.GlobalEnv$dacPos]+numberDAC
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
colnames(res1)=c("time",y_names)

ggplot(res1,aes(x=time))+geom_line(aes(y=EBV))
ggplot(res1,aes(x=time))+geom_line(aes(y=NK))
ggplot(res1,aes(x=time))+geom_line(aes(y=ODC_le1))


