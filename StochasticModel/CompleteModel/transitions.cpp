// Parameters list:

double TeE = 0.4 ;
double TrE = 0.2 ;
double Tr2 = 0.09 ;  
double Te2 = 0.5 ;

double cIL2 = 200;
double cMem = 200;
double cEBV = 1000;
double cDAC= 13.02883;
double nk2 =  0.04166667; //....-> 1/24
double probDup= 0.6666667;//....-> 2/3

  
static map <string, double> InjEBVTime={{"FirstInj", 168 },{"SecondInj", 1608 },{"ThirdInj", 3048 },{"FourInj", 4488 },{"FifthInj", 5928 }}; // hour scale
static map <string, double> RATES_killingTransitions={{"TeffkillsEBV", 0.1 },{"TregKillsTeff", 1 },{"NKkillsTreg", 0.1 },{"NKkillsTeff", 0.1 },{"TeffKillsODC", 0.15  }};

// General transitions:

//////////////////////////////////////////////////////////////////////////////////////
// transitions from Timo to Tcells and NK
//////////////////////////////////////////////////////////////////////////////////////

double TimoReg(double *Value, map <string,int>& NumTrans, map <string,int>& NumPlaces,const vector<string> & NameTrans, const struct InfTr* Trans, const int T, const double& time)
{
    double rate=0.0 ;
    int idx = NumPlaces.find("Resting_Treg_temp") -> second ;
    
    double RestTregOut=0.0 ;
    RestTregOut = Value[idx];
    
    rate = ( RestTregOut / 63.0 ) * 20.0;

	return rate;
}

double TimoEff(double *Value, map <string,int>& NumTrans, map <string,int>& NumPlaces,const vector<string> & NameTrans, const struct InfTr* Trans, const int T, const double& time)
{
    double rate=0.0;
    int idx = NumPlaces.find("Resting_Teff_temp") -> second ;
    
    double RestTeffOut=0.0 ;
    RestTeffOut = Value[idx];
    
    rate = ( RestTeffOut/1687.0 )* 500.0;
     
    return rate;
}

double NKBorn(double *Value, map <string,int>& NumTrans, map <string,int>& NumPlaces,const vector<string> & NameTrans, const struct InfTr* Trans, const int T, const double& time)
{
    double rate=0.0;
    int idx = NumPlaces.find("NK_temp") -> second ;
    
    double NKOut=0.0 ;
    NKOut = Value[idx];    
    
    rate = (NKOut/375 )* 100.0;
    
	return rate;
}

//////////////////////////////////////////////////////////////////////////////////////
// Duplication of the Tcells
//////////////////////////////////////////////////////////////////////////////////////

double TeffDup(double *Value, map <string,int>& NumTrans, map <string,int>& NumPlaces,const vector<string> & NameTrans, const struct InfTr* Trans, const int T, const double& time)
{
    double rate=0.0;
    int idx = NumPlaces.find("DAC") -> second ;
    int idx2 = NumPlaces.find("IL2") -> second ;
    int idx3 = NumPlaces.find("Teff") -> second ;
   
    double DAC=0.0 ;
    double IL2=0.0 ;
    double Teff=0.0 ;
    double CellTot=1 ;
    
    DAC = Value[idx];       
    IL2 = Value[idx2];        
    Teff = Value[idx3];
    
    double sum =0.0;
	for (auto it=NumPlaces.begin(); it!=NumPlaces.end(); it++)
    {
        if( it-> first != "EffectorMemory" && it-> first != "NK_temp" && it-> first != "Resting_Teff_temp"&& it-> first != "Resting_Treg_temp")
            sum  += Value[it-> second];
    }
    
    CellTot = sum;
    
    rate = Te2 * ( 1-exp(-IL2/cIL2)  )* ( exp(-DAC/cDAC) ) * Teff * IL2 / CellTot ;
    
    if(NameTrans[T]=="TeffDup_Asym")
        return (1-probDup)*rate;
    else return probDup*rate;
    
}

double TregDup(double *Value, map <string,int>& NumTrans, map <string,int>& NumPlaces,const vector<string> & NameTrans, const struct InfTr* Trans, const int T, const double& time)
{
    double rate=0.0;
    int idx = NumPlaces.find("DAC") -> second ;
    int idx2 = NumPlaces.find("IL2") -> second ;
    int idx3 = NumPlaces.find("Treg") -> second ;
   
    double DAC=0.0 ;
    double IL2=0.0 ;
    double Treg=0.0 ;
    double CellTot=1.0 ;
    
    DAC = Value[idx];       
    IL2 = Value[idx2];        
    Treg = Value[idx3];
    
    double sum =0.0;
	for (auto it=NumPlaces.begin(); it!=NumPlaces.end(); it++)
    {
        if( it-> first != "EffectorMemory" && it-> first != "NK_temp" && it-> first != "Resting_Teff_temp"&& it-> first != "Resting_Treg_temp")
            sum  += Value[it-> second];
    }
    
    CellTot = sum;
    
    rate = Tr2 * ( 1.0 -exp(-IL2/cIL2)  )* ( exp(-DAC/cDAC) ) * Treg * IL2 / CellTot ;

    return rate;
    
}

double NKdup(double *Value, map <string,int>& NumTrans, map <string,int>& NumPlaces,const vector<string> & NameTrans, const struct InfTr* Trans, const int T, const double& time)
{
    double rate=0.0;
    int idx = NumPlaces.find("NK") -> second ;
    int idx2 = NumPlaces.find("IL2") -> second ;
   
    double NK=0.0 ;
    double IL2=0.0 ;
    double CellTot=1.0 ;
    
    NK = Value[idx];   
    IL2 = Value[idx2];  
    
    double sum =0.0;
	for (auto it=NumPlaces.begin(); it!=NumPlaces.end(); it++)
    {
        string s;
        s = it-> first;
        if( s != "EffectorMemory" && s != "NK_temp" && s != "Resting_Teff_temp" && s != "Resting_Treg_temp")
        {
            sum  += Value[it-> second];
        }
    }
    
    CellTot = sum;
    
    rate = nk2 * NK * IL2 / CellTot ;

    
    return rate;
    
}

//////////////////////////////////////////////////////////////////////////////////////
// Activation of the Tcells
//////////////////////////////////////////////////////////////////////////////////////

double TeffActivation(double *Value, map <string,int>& NumTrans, map <string,int>& NumPlaces,const vector<string> & NameTrans, const struct InfTr* Trans, const int T, const double& time)
{

double rate=0.0;
    int idx = NumPlaces.find("EBV") -> second ;
    int idx2 = NumPlaces.find("Resting_Teff") -> second ;
   
    double EBV=0.0 ;
    double Resting_Teff=0.0;
     
    EBV = Value[idx];       
    Resting_Teff = Value[idx2]; 
    
    rate = TeE * ( 1.0 - exp(-EBV/cEBV)  ) * Resting_Teff;  

    return rate;
    
}

double TregActivation(double *Value, map <string,int>& NumTrans, map <string,int>& NumPlaces,const vector<string> & NameTrans, const struct InfTr* Trans, const int T, const double& time)
{    
    double rate=0.0;
    int idx = NumPlaces.find("EBV") -> second ;
    int idx2 = NumPlaces.find("Resting_Treg") -> second ;
    int idx3 = NumPlaces.find("Teff") -> second ;
   
    double EBV=0.0 ;
    double Teff=0.0 ;
    double Resting_Treg=0.0;
     
    EBV = Value[idx];       
    Resting_Treg = Value[idx2]; 
    Teff = Value[idx3];  
    
    rate = TrE * ( (Teff )/(Teff + EBV + 1.0) ) * Resting_Treg;

    return rate;

}


double MemActivation(double *Value, map <string,int>& NumTrans, map <string,int>& NumPlaces,const vector<string> & NameTrans, const struct InfTr* Trans, const int T, const double& time)
{
    double MemE = 0.0;
    double rate=0.0;
    int idx = NumPlaces.find("EBV") -> second ;
    int idx2 = NumPlaces.find("EffectorMemory") -> second ;
    
    double EBV=0.0 ;
    double Mem=0.0 ;
    
    EBV = Value[idx];       
    Mem = Value[idx2]; 
    
    ifstream in;    // Create an input file stream.
    in.open("SecondInj.txt");  // Use it to read from a file named data.txt.
    int x;
    in >> x;   
        
    if(x == 1) 
        MemE = 2.0 * TeE * ( 1.0 -exp(-Mem/cMem)  ) ;
  
    rate = MemE * ( 1.0 - exp(-EBV/cEBV) ) * Mem;
        
    return rate;
}

//////////////////////////////////////////////////////////////////////////////////////
// Regulation of the Tcells
//////////////////////////////////////////////////////////////////////////////////////

double Killing(double *Value, map <string,int>& NumTrans, map <string,int>& NumPlaces,const vector<string> & NameTrans, const struct InfTr* Trans, const int T, const double& time)
{

    double rate=0.0;
    double CellTot=1.0;
    double sum =0.0;
    
	for (auto it=NumPlaces.begin(); it!=NumPlaces.end(); it++)
    {
        if( it-> first != "EffectorMemory" && it-> first != "NK_temp" && it-> first != "Resting_Teff_temp"&& it-> first != "Resting_Treg_temp")
            sum  += Value[it-> second];
    }
    
    CellTot = sum;
    
    rate = RATES_killingTransitions.find(NameTrans[T]) -> second ;
    
    double intensity = 1.0;
    for (unsigned int k=0; k<Trans[T].InPlaces.size(); k++)
    {
        intensity *= pow(Value[Trans[T].InPlaces[k].Id],Trans[T].InPlaces[k].Card);
    }
    
    return rate*intensity / CellTot;    
}


double KillingODC(double *Value, map <string,int>& NumTrans, map <string,int>& NumPlaces,const vector<string> & NameTrans, const struct InfTr* Trans, const int T, const double& time)
{

    double rate=0.0;
    double CellTot=1.0;
    double sum =0.0;
    
	for (auto it=NumPlaces.begin(); it!=NumPlaces.end(); it++)
    {
        if( it-> first != "EffectorMemory" && it-> first != "NK_temp" && it-> first != "Resting_Teff_temp"&& it-> first != "Resting_Treg_temp")
            sum  += Value[it-> second];
    }
    
    CellTot = sum;

    rate = RATES_killingTransitions.find("TeffKillsODC") -> second ;
    
    double intensity = 1.0;
    for (unsigned int k=0; k<Trans[T].InPlaces.size(); k++)
    {
        intensity *= pow(Value[Trans[T].InPlaces[k].Id],Trans[T].InPlaces[k].Card);
    }
    
    return rate*intensity / CellTot;    
}
