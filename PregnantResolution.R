####################### Resolution Pregnant Model:

library(ggplot2)
library(plyr)

############################################
#### Inj times

numberEBV <<- 1000

InjEBVTime=c(7,60,110,210,310,410,510)*24
InjPregTime= c(80,170,260,350)*24


L<-length(c(InjPregTime,InjEBVTime))

InjTime<-data.frame("Times" = numeric(L), "kind" =character(L) )
InjTime$Times = c(InjPregTime,InjEBVTime)
InjTime$kind=rep(c("Preg","EBV"),c(length(InjPregTime),length(InjEBVTime)) )
InjTime<-InjTime[with(InjTime, order(Times)),]
rep<-count(InjTime[,1])

if( TRUE%in%c(rep[,2]>1) )
{
  rep_pos<-which(rep[rep[,2]>1,1]==InjTime)
  InjTime$kind[rep_pos]="BOTH"
  InjTime<-unique(InjTime)
}
############################################

FinalTime<-730*24

step=24

## Parameters for sick woman
# Treg Activation from .2 -> .8
# Teff Activation from .4 -> .1

s<-runif(100,min = 0.05,.11) # generation of the 100 trajectories


for (j in 1:length(s)) {
  preg.index<-0

p0<-c(.4,.2,0.09,0.5,0.15, 1 ,0.1 ,0.1,0.1)
p1<-c(.4- s[j],.2+ 2*s[j],0.09,0.5,0.15, 1 ,0.1 ,0.1,0.1)
p2<-c(.4- 2*s[j],.2+ 4* s[j],0.09,0.5,0.15, 1 ,0.1 ,0.1,0.1)
p3<-c(.4- 3*s[j],.2+ 6*s[j],0.09,0.5,0.15, 1 ,0.1 ,0.1,0.1)

p.all<-data.frame(p0=p0,p1=p1,p2=p2,p3=p3)

row.names(p.all)<- c( "TeE","TrE","Tr2","Te2","TekODC","TrkTe","TekEBV","rec","NKkT")


source("ModelRRMS.R")

p<-p0
names(p)<- c( "TeE","TrE","Tr2","Te2","TekODC","TrkTe","TekEBV","rec","NKkT")

resPreg <-lsoda(yini,seq(from = 0, to = InjTime[1,1], by = step), funODE, parms=p)

for(i in 1:length(InjTime[,1]))
{
  yini<-tail(resPreg,1)[,-1]
  
  
  
  if(InjTime[i,2]=="EBV"){
    ebvPosition<-sample(.GlobalEnv$ebvPos,1)
    yini[ebvPosition]<-yini[ebvPosition]+numberEBV
  }else  if(InjTime[i,2]=="Preg"){
    preg.index=preg.index+1
  }else{
    preg.index=preg.index+1
    ebvPosition<-sample(.GlobalEnv$ebvPos,1)
    yini[ebvPosition]<-yini[ebvPosition]+numberEBV
  }
  
  #### check the pregnancy
  if(preg.index!=0 && preg.index!=4) {
    if(preg.index==1) p <- p.all[,2] 
    if(preg.index==2) p <- p.all[,3] 
    if(preg.index==3) p <- p.all[,4] 
  }else p <- p.all[,1]
    
    names(p)<- c( "TeE","TrE","Tr2","Te2","TekODC","TrkTe","TekEBV","rec","NKkT")
  
  start_time<-InjTime[i,1]
  
  if(i==length(InjTime[,1]))
  {
    end_time<-FinalTime
  }else{
    end_time<-InjTime[(i+1),1]
  }
  #res1[length(res1[,2]),1]<- tail(res1[,1],1)-.2
  resPreg <- resPreg[-length(resPreg[,1]),]
  res <-lsoda(yini,seq(from = start_time, to = end_time, by = step),funODE,parms=p)
  
  resPreg<-rbind(resPreg,res)

}

resPreg<-cbind(resPreg,proportion=s[j])
write.table(file=paste(j,"ODE.txt",sep=""),resPreg)

}



####################################
## PLOTS


filesAll <- list.files('./',pattern='^[0-9]{1,4}[A-Z]{3}.txt')
result<-list()

for(i in 1:length(filesAll) ){ 
  result[[i]] <-read.csv(paste('./',filesAll[i],sep=""), sep="")
}

dati<-ldply(result, data.frame)

dati<- as.data.frame(resPreg)

ggplot(dati,aes(x=time/24))+geom_line(aes(y=ODC_le1,group=proportion,colour=proportion))

ggplot(dati,aes(x=time/24))+geom_line(aes(y=ODC_le1_px1_py1_pz1,group=proportion,colour=proportion))+theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"),legend.text=element_text(size=18),legend.title=element_text(size=20,face="bold")) + ylab("ODC cells \nirreversibly damaged") + xlab("Days")+ labs(col = "Treg-Teff \nbalance")+ scale_y_continuous(sec.axis = dup_axis(trans = ~./500, name = "Death percetage",labels = scales::percent  ) ) + geom_vline(xintercept=c(80,170,260,350),col="black") + geom_line(data=as.data.frame(resPreg),aes(x=time/24,y=ODC_le1_px1_py1_pz1),col="red")+xlim(0,500)

ggplot(dati,aes(x=time/24))+geom_line(aes(y=EBV_px1_py1_pz1,group=proportion,colour=proportion))+theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"),legend.text=element_text(size=18),legend.title=element_text(size=20,face="bold")) + ylab("ODC cells \nirreversibly damaged") + xlab("Days")+ labs(col = "Treg-Teff shift")+ geom_vline(xintercept=c(80,170,260,350),col="red") + geom_line(data=as.data.frame(resPreg),aes(x=time/24,y=EBV_px1_py1_pz1),col="red")

ggplot(dati,aes(x=time/24))+geom_line(aes(y=Teff_px1_py1_pz1,group=proportion,colour=proportion))+theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"),legend.text=element_text(size=18),legend.title=element_text(size=20,face="bold")) + ylab("Teff cells") + xlab("Days")+ labs(col = "Treg-Teff shift") + geom_vline(xintercept=c(80,170,260,350),col="red") + geom_line(data=as.data.frame(resPreg),aes(x=time/24,y=Teff_px1_py1_pz1),col="red")


ggplot(dati,aes(x=time/24))+geom_line(aes(y=NK_px1_py1_pz1,group=proportion,colour=proportion))+theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"),legend.text=element_text(size=18),legend.title=element_text(size=20,face="bold")) + ylab("Teff cells") + xlab("Days")+ labs(col = "Treg-Teff shift") + geom_vline(xintercept=c(80,170,260,350),col="red") + geom_line(data=as.data.frame(resPreg),aes(x=time/24,y=NK_px1_py1_pz1),col="red")

###########################

resPreg1<-dati[dati$proportion==max(dati$proportion),]


ggplot(as.data.frame(resPreg),aes(x=time/24))+geom_line(aes(y=Teff_px1_py1_pz1,colour="Not Pregnant"))+ geom_line(data=resPreg1,aes(x=time/24,y=Teff_px1_py1_pz1,colour="Pregnant"))+ theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"),legend.text=element_text(size=18),legend.title=element_text(size=20,face="bold")) + ylab("Teff cells") + xlab("Days")+ labs(col = "Woman with RRMS") + geom_vline(xintercept=c(80,170,260,350),col="black") 

ggplot(as.data.frame(resPreg),aes(x=time/24))+geom_line(aes(y=Treg_px1_py1_pz1,colour="Not Pregnant"))+ geom_line(data=resPreg1,aes(x=time/24,y=Treg_px1_py1_pz1,colour="Pregnant"))+ theme(axis.text=element_text(size=18),axis.title=element_text(size=20,face="bold"),legend.text=element_text(size=18),legend.title=element_text(size=20,face="bold")) + ylab("Treg cells") + xlab("Days")+ labs(col = "Woman with RRMS") + geom_vline(xintercept=c(80,170,260,350),col="black") 

ggplot(resPreg,aes(x=time/24))+geom_line(aes(y=resPreg[,teffPos+1],col="Teff") )+geom_line(aes(y=resPreg[,tregPos+1],col="Treg"))

ggplot(resPreg,aes(x=time))+geom_line(aes(y=ODC_le1_px1_py1_pz1))+geom_vline(xintercept=c(80,170,260,350),col="red")


