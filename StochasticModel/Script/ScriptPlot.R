setwd("~/Desktop/CMISF2019/NOArcs/")
library("sm")
library(ggplot2)
library(readr)

# 100 trajectories
#
# Results are saved in: TauSolutions
# 
# =========================== TIME ===========================
#   
#   Total time required: 0s.
# 
# 
# Results are saved in: SSAsolutions
# 
# =========================== TIME ===========================
#   
#   Total time required: 8s


SSAsolutions <- read_table2("SSASolutions.trace")
TauSolutions <- read_table2("TAUSolutions.trace")

unique(table(SSAsolutions$Time))->lengthSSA
unique(table(TauSolutions$Time))->lengthTAU

SSAsolutions$ID<-rep(1:lengthSSA,each=length(table(SSAsolutions$Time)))
TauSolutions$ID<-rep(1:lengthTAU,each=length(table(TauSolutions$Time)))

ggplot(SSAsolutions,aes(x=Time,y=X1,group=ID))+geom_line()+labs(title="SSA")
ggplot(TauSolutions,aes(x=Time,y=X1,group=ID))+geom_line()+labs(title="TAU")



########################################################################
########### Histogram difference

SSASol <- unlist(read_csv("SSASol.mtx"))
TAUSol_15 <- unlist(read_csv("TAUSol_15.mtx"))
TAUSol_2 <- unlist(read_csv("TAUSol_2.mtx"))
TAUSol_25 <- unlist(read_csv("TAUSol_25.mtx"))
TAUSol_3 <- unlist(read_csv("TAUSol_3.mtx"))
TAUSol_07<-unlist(read_csv("TAUSol_07.mtx"))
TAUSol_06<-unlist(read_csv("TAUSol_06.mtx"))
TAUSol_05<-unlist(read_csv("TAUSol_05.mtx"))
TAUSol_04<-unlist(read_csv("TAUSol_04.mtx"))

points<-seq(0,1000,2)

a25<-sm.density(TAUSol_25,eval.points=points)
a15<-sm.density(TAUSol_15,eval.points=points)
a2<-sm.density(TAUSol_2,eval.points=points)
a3<-sm.density(TAUSol_3,eval.points=points)
a05<-sm.density(TAUSol_05,eval.points=points)
a06<-sm.density(TAUSol_06,eval.points=points)
a07<-sm.density(TAUSol_07,eval.points=points)
a04<-sm.density(TAUSol_04,eval.points=points)

b<-sm.density(SSASol,eval.points=points)


err_25<-sum(abs(a25$estimate-b$estimate))/2
err_15<-sum(abs(a15$estimate-b$estimate))/2
err_2<-sum(abs(a2$estimate-b$estimate))/2
err_3<-sum(abs(a3$estimate-b$estimate))/2
err_04<-sum(abs(a04$estimate-b$estimate))/2
err_05<-sum(abs(a05$estimate-b$estimate))/2
err_06<-sum(abs(a06$estimate-b$estimate))/2
err_07<-sum(abs(a07$estimate-b$estimate))/2

err<-c(err_04,err_05,err_06,err_07,err_15,err_2,err_25,err_3)
eps<-c(.04,.05,.06,.07,.15,.2,.25,.3)
  
err.df<-data.frame(eps=eps,err=err)


ggplot(err.df, aes(x=eps,y=err)) + geom_line()+geom_point()
