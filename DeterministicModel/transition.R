## Parameters list:

TeE = 0.4 ;
TrE = 0.2 ;
Tr2 = 0.09 ;  
Te2 = 0.5 ;

cIL2 = 200;
cMem = 200;
cEBV = 1000;
cDAC= 13.02883;
nk2 =  1/24
probDup=  2/3

InjEBVTime=  c(168 ,1608, 3048, 4488, 5928)  # hour scale
names(InjEBVTime)=c("FirstInj", "SecondInj","ThirdInj", "FourthInj","FifthInj")

RATES_killingTransitions=c(0.1 ,1,0.1, 0.1 , 0.15)
names(RATES_killingTransitions)=c( "TeffkillsEBV","TregKillsTeff","NKkillsTreg","NKkillsTeff" ,"TeffKillsODC") 


# General transitions:

###########################################
# transitions from Timo to Tcells, and NK
###########################################

TimoReg<-function(Value,NameTrans, NumPlaces, time,InputPlaces,card)
{
    rate=0.0 ;
    idx = NumPlaces["Resting_Treg_temp"]
    
    RestTregOut=0.0 ;
    RestTregOut = Value[idx];
    
    rate = ( RestTregOut / 63.0 ) * 20.0;

	return(rate)
}

 TimoEff<-function(Value,NameTrans, NumPlaces, time,InputPlaces,card)
{
     rate=0.0;
    idx = NumPlaces["Resting_Teff_temp"] ;
    
     RestTeffOut=0.0 ;
    RestTeffOut = Value[idx];
    
    rate = ( RestTeffOut/1687.0 )* 500.0;
     
    return(rate)
}

 NKBorn<-function(Value,NameTrans, NumPlaces, time,InputPlaces,card)
{
     rate=0.0;
    idx = NumPlaces["NK_temp"] ;
    
     NKOut=0.0 ;
    NKOut = Value[idx];    
    
    rate = (NKOut/375 )* 100.0;
    
	return(rate)
}

###########################################
# Duplication of the Tcells
###########################################

 TeffDup<-function(Value,NameTrans, NumPlaces, time,InputPlaces,card)
{
     rate=0.0;
    idx = NumPlaces["DAC"] ;
    idx2 = NumPlaces["IL2"] ;
    idx3 = NumPlaces["Teff"] ;
   
     DAC=0.0 ;
     IL2=0.0 ;
     Teff=0.0 ;
     CellTot=1 ;
    
    DAC = Value[idx];       
    IL2 = Value[idx2];        
    Teff = Value[idx3];
    
     sum =0.0;
	for (it in NumPlaces)
    {
        if( names(NumPlaces[it]) != "EffectorMemory" &&  names(NumPlaces[it]) != "NK_temp" && names(NumPlaces[it]) !="Resting_Teff_temp"&& names(NumPlaces[it]) !="Resting_Treg_temp")
            sum = sum + Value[it];
    }
    
    CellTot = sum;
    
    rate = Te2 * ( 1-exp(-IL2/cIL2)  )* ( exp(-DAC/cDAC) ) * Teff * IL2 / CellTot ;
    
    if(NameTrans=="TeffDup_Asym")
        return( (1-probDup)*rate)
    else return(probDup*rate)
    
}

 TregDup<-function(Value,NameTrans, NumPlaces, time,InputPlaces,card)
{
     rate=0.0;
    idx = NumPlaces["DAC"] ;
    idx2 = NumPlaces["IL2"] ;
    idx3 = NumPlaces["Treg"] ;
   
     DAC=0.0 ;
     IL2=0.0 ;
     Treg=0.0 ;
     CellTot=1.0 ;
    
    DAC = Value[idx];       
    IL2 = Value[idx2];        
    Treg = Value[idx3];
    
     sum =0.0;
  for (it in NumPlaces)
  {
    if( names(NumPlaces[it]) != "EffectorMemory" &&  names(NumPlaces[it]) != "NK_temp" && names(NumPlaces[it]) !="Resting_Teff_temp"&& names(NumPlaces[it]) !="Resting_Treg_temp")
      sum = sum + Value[it];
  }
    
    CellTot = sum;
    
    rate = Tr2 * ( 1.0 -exp(-IL2/cIL2)  )* ( exp(-DAC/cDAC) ) * Treg * IL2 / CellTot ;

    return(rate)
    
}

 NKdup<-function(Value,NameTrans, NumPlaces, time,InputPlaces,card)
{
     rate=0.0;
    idx = NumPlaces["NK"] ;
    idx2 = NumPlaces["IL2"] ;
   
     NK=0.0 ;
     IL2=0.0 ;
     CellTot=1.0 ;
    
    NK = Value[idx];   
    IL2 = Value[idx2];  
    
     sum =0.0;
  for (it in NumPlaces)
  {
    if( names(NumPlaces[it]) != "EffectorMemory" &&  names(NumPlaces[it]) != "NK_temp" && names(NumPlaces[it]) !="Resting_Teff_temp"&& names(NumPlaces[it]) !="Resting_Treg_temp")
      sum = sum + Value[it];
  }
    
    CellTot = sum;
    
    rate = nk2 * NK * IL2 / CellTot ;

    
    return(rate)
    
}

###########################################
# Activation of the Tcells
###########################################

 TeffActivation<-function(Value,NameTrans, NumPlaces, time,InputPlaces,card)
{

 rate=0.0;
    idx = NumPlaces["EBV"] ;
    idx2 = NumPlaces["Resting_Teff"] ;
   
     EBV=0.0 ;
     Resting_Teff=0.0;
     
    EBV = Value[idx];       
    Resting_Teff = Value[idx2]; 
    
    rate = TeE * ( 1.0 - exp(-EBV/cEBV)  ) * Resting_Teff;  

    return(rate)
    
}

 TregActivation<-function(Value,NameTrans, NumPlaces, time,InputPlaces,card)
{    
     rate=0.0;
    idx = NumPlaces["EBV"] ;
    idx2 = NumPlaces["Resting_Treg"] ;
    idx3 = NumPlaces["Teff"] ;
   
     EBV=0.0 ;
     Teff=0.0 ;
     Resting_Treg=0.0;
     
    EBV = Value[idx];       
    Resting_Treg = Value[idx2]; 
    Teff = Value[idx3];  
    
    rate = TrE * ( (Teff )/(Teff + EBV + 1.0) ) * Resting_Treg;

    return(rate)

}


 MemActivation<-function(Value,NameTrans, NumPlaces, time,InputPlaces,card)
{
     MemE = 0.0;
     rate=0.0;
    idx = NumPlaces["EBV"] ;
    idx2 = NumPlaces["EffectorMemory"] ;
    
     EBV=0.0 ;
     Mem=0.0 ;
    
    EBV = Value[idx];       
    Mem = Value[idx2]; 
    
    if(time > as.integer(InjEBVTime["SecondInj"] )) 
      MemE = 2.0 * TeE * ( 1.0 -exp(-Mem/cMem)  ) ;
  
    rate = MemE * ( 1.0 - exp(-EBV/cEBV) ) * Mem;
        
    return(rate)
}

###########################################
# Regulation of the Tcells
###########################################

 Killing<-function(Value,NameTrans, NumPlaces, time,InputPlaces,card)
{

     rate=0.0;
     CellTot=1.0;
     
     sum =0.0;
     for (it in NumPlaces)
     {
       if( names(NumPlaces[it]) != "EffectorMemory" &&  names(NumPlaces[it]) != "NK_temp" && names(NumPlaces[it]) !="Resting_Teff_temp"&& names(NumPlaces[it]) !="Resting_Treg_temp")
         sum = sum + Value[it];
     }
    
    CellTot = sum;
    
    rate = RATES_killingTransitions[paste(NameTrans)]
    
    intensity = 1.0;
    for ( k in 1:length(InputPlaces))
    {
      indexx<-NumPlaces[InputPlaces[k]]
      intensity = intensity * Value[indexx]^(card[k]);
    }  
    
    return(rate*intensity / CellTot)    
}


 KillingODC<-function(Value,NameTrans, NumPlaces, time,InputPlaces,card)
{

     rate=0.0;
     CellTot=1.0;

     sum =0.0;
     for (it in NumPlaces)
     {
       if( names(NumPlaces[it]) != "EffectorMemory" &&  names(NumPlaces[it]) != "NK_temp" && names(NumPlaces[it]) !="Resting_Teff_temp"&& names(NumPlaces[it]) !="Resting_Treg_temp")
         sum = sum + Value[it];
     }
    CellTot = sum;

    rate = RATES_killingTransitions["TeffKillsODC"] ;
    
     intensity = 1.0;
    for ( k in 1:length(InputPlaces))
    {
      indexx<-NumPlaces[InputPlaces[k]]
        intensity = intensity * Value[indexx]^(card[k]);
    }
    
    return( rate*intensity / CellTot)    
}
