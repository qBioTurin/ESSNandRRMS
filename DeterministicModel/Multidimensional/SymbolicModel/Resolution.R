library(deSolve)
library(ggplot2)
library(plyr)
setwd("~/Desktop/CIBBproceeding/Simmetrie/")

############################################
#### Inj times ########
#DACther=T
InjEBVTime=c(7,67,127,187,247)*24

if(DACther)
{
  InjDACTime= seq(67,337,by=30)*24
  
}else InjDACTime= NULL

numberDAC <<- 300
numberEBV <<- 10000


EBVplaces.injection<-"P1"
DACplaces.injection<-"P2"

L<-length(c(InjDACTime,InjEBVTime))

InjTime<-data.frame("Times" = numeric(L), "kind" =character(L) )
InjTime$Times = c(InjDACTime,InjEBVTime)
InjTime$kind=rep(c("DAC","EBV"),c(length(InjDACTime),length(InjEBVTime)) )
InjTime<-InjTime[with(InjTime, order(Times)),]
rep<-count(InjTime[,1])

if( TRUE%in%c(rep[,2]>1) )
{
	rep_pos<-which(InjTime$Times%in%rep[rep[,2]>1,1])
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

source("Griglia3dNEW.R")

# yini["EBV_P1"]=10000
# resSimm <-lsoda(yini,0:11, funODE, parms=p)

resSimm <-lsoda(yini,seq(from = 0, to = InjTime[1,1] , by = step), funODE, parms=p)


for(i in 1:length(InjTime[,1]))
{
	yini<-tail(resSimm,1)[,-1]

	if(InjTime[i,2]=="EBV"){
		yini["EBV_P1"]<-yini["EBV_P1"]+numberEBV
	}else  if(InjTime[i,2]=="DAC"){
		yini["DAC_P2"]<-yini["DAC_P2"]+numberDAC
	}else{
		yini["DAC_P2"]<-yini["DAC_P2"]+numberDAC
		yini["EBV_P1"]<-yini["EBV_P1"]+numberEBV
	}

	start_time<-InjTime[i,1]

	if(i==length(InjTime[,1]))
	{
		end_time<-FinalTime
	}else{
		end_time<-InjTime[(i+1),1]
	}
	#res1[length(res1[,2]),1]<- tail(res1[,1],1)-.2
	resSimm <- resSimm[-length(resSimm[,1]),]
	
	#res <-lsoda(yini,168:180,funODE,parms=p)
	
	res <-lsoda(yini,seq(from = start_time, to = end_time, by = step),funODE,parms=p)

	resSimm<-rbind(resSimm,res)
}

resSimm<-as.data.frame(resSimm)
resSimm[,"time"]<-resSimm[,"time"]/24




