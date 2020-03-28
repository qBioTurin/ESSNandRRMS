
moltColorDom = c(2,2,23)
names(moltColorDom) =c("P1","P2","P3")

# costant definition
cIL2 <<- 200
cMem <<- 200
cEBV <<- 1000
cTcell <- 200
cTeff <- 200
#cDAC <<- - numberDAC /log(.1)
cDAC <<- - 1 /log(.1)

probDup<<- 2/3

## injections at time fixed
DACinj = 0 ## it is not MA
EBVs = 0  ## it is not MA


NKdup=1/24
TregDeath=1/24
TeffDeath=1/24
NKDegradation=1/24

NKarrive=100/375
FromTimoREG=20/63
FromTimoEFF=500/1687
MovementEBV=.1
DACMovements=.1
MovementTeff=.1
MovementTreg=.1

funODE <- function(t,y, parms)
{
	y[y<1e-8]=0
	names(y)=y_names

	## Parameter to estimate
	TeE = parms["TeE"]
	TregActivation = parms["TrE"]
	TregDup = parms["Tr2"]
	Te2 = parms["Te2"]
	TekODC = parms["TekODC"]
	TregKillsTeff = parms["TrkTe"]
	TeffkillsEBV = parms["TekEBV"]
	Remyelinization = parms["rec"]
	NKkillT=parms["NKkT"]

	if( is.na(parms["DACD"]) )
	{
		DACDegradation = 0
		cDAC<<-1
	}
	else DACDegradation = parms["DACD"]

	### GENERAL TRANSITIONS:

	mP1=moltColorDom["P1"]
	mP2=moltColorDom["P2"]
	mP3=moltColorDom["P3"]

	EBVgrid<-y[ebvPos]
	Teffgrid<-y[teffPos]
	Treggrid<-y[tregPos]

	sumEBVgrid <- y["EBV_P1"]*mP1 + y["EBV_P2"]*mP2 +y["EBV_P3"]*mP3
	sumTeffgrid  <- y["Teff_P1"]*mP1 + y["Teff_P2"]*mP2 +y["Teff_P3"]*mP3
	sumTcellgrid <- sum(y[c("Teff_P1","Teff_P1")]*mP1 + y[c("Teff_P2","Teff_P2")]*mP2 +y[c("Teff_P3","Teff_P3")]*mP3)

	if(sumEBVgrid == 0 ){ p_eff <- rep(0,length(EBVgrid)) }else{ p_eff<- EBVgrid/sumEBVgrid }

	if(sumTeffgrid == 0 ){ p_reg <- rep(0,length(Teffgrid)) }else{ p_reg<- Teffgrid/sumTeffgrid }

	if(sumTcellgrid == 0 ){ p_dac <- rep(0,length(Teffgrid)) }else{ p_dac <- (Treggrid+Teffgrid)/sumTcellgrid }

	names(p_eff)<- combos
	names(p_reg)<- combos
	names(p_dac)<- combos

	p_ebv<-1/(mP1+mP2+mP3)

	for (position in combos )
	{
		m <- moltColorDom[position]

		IL2 <-y[paste0("IL2_",position)]
		Treg<-y[paste0("Treg_",position)]
		Teff<-y[paste0("Teff_",position)]
		EBV<-y[paste0("EBV_",position)]
		NK<-y[paste0("NK_",position)]
		DAC<-y[paste0("DAC_",position)]

		ODCl1<-y[paste0("ODC_L1_",position)]
		ODCl2<-y[paste0("ODC_L2_",position)]
		ODCl3<-y[paste0("ODC_L3_",position)]
		ODCl4<-y[paste0("ODC_L4_",position)]
		ODCl5<-y[paste0("ODC_L5_",position)]

		ODC<-ODCl1+ODCl2+ODCl3+ODCl4+ODCl5

		Mem<-y[memPos]
		RestingTreg<-y[paste0("Resting_Treg_",position)]
		RestingTeff<-y[paste0("Resting_Teff_",position)]

		TotMeet_xyz<- ODC+Teff+Treg+DAC+IL2+NK+RestingTreg+RestingTeff+EBV

		####### parameters marking depending calculation

		if(t <= InjEBVTime[2] ){MemE=0}else{ MemE = 2 * TeE * ( 1-exp(-Mem/cMem)  )}

		###############################################

		R_TeffDup<- 1/TotMeet_xyz*  ( 1-exp(-IL2/cIL2)  )* ( exp(-DAC/cDAC) ) * IL2 * Teff
		assign(paste("R_TeffDup_Sym_p1", position, sep="_"),  probDup * R_TeffDup )
		assign(paste("R_TeffDup_Asym_p1", position, sep="_"),  (1-probDup) * R_TeffDup )
		assign(paste("R_TregDup_p1", position, sep="_"),1/TotMeet_xyz*  ( 1-exp(-IL2/cIL2)  )* ( exp(-DAC/cDAC) ) * IL2 * Treg )
		assign(paste("R_NKdup_p1", position, sep="_"),1/TotMeet_xyz* ( exp(-DAC/cDAC) ) * IL2 * NK )

		assign(paste("R_TeffActivation_p1", position, sep="_"),  ( 1 - exp(-EBV/cEBV))* RestingTeff )
		assign(paste("R_TregActivation_p1", position, sep="_"),  ( (Teff )/(Teff + EBV+ 1) ) * RestingTreg )
		assign(paste("R_MemActivation_p1", position, sep="_"),  ( 1- exp(-EBV/cEBV) ) * Mem )

		assign(paste("R_TeffkillsEBV_p1", position, sep="_"),1/TotMeet_xyz*EBV*Teff )
		assign(paste("R_TregKillsTeff_p1", position, sep="_"),1/TotMeet_xyz*Treg*Teff)
		assign(paste("R_NKkillsTreg_p1", position, sep="_"),1/TotMeet_xyz*Treg*NK )
		assign(paste("R_NKkillsTeff_p1", position, sep="_"),1/TotMeet_xyz*Teff*NK )

		assign(paste("R_TeffKillsODC1_l1_L1_p1", position, sep="_"), 1/TotMeet_xyz*ODCl2*Teff )
		assign(paste("R_TeffKillsODC2_l1_L2_p1", position, sep="_"), 1/TotMeet_xyz*ODCl3*Teff )
		assign(paste("R_TeffKillsODC3_l1_L3_p1", position, sep="_"), 1/TotMeet_xyz*ODCl4*Teff )
		assign(paste("R_TeffKillsODC4_l1_L4_p1", position, sep="_"), 1/TotMeet_xyz*ODCl5*Teff )

		# Movements:
		# for istance: R_DACMovements_p1_P1_p2_P1

		for( indexpos in 1:length(combos))
		{
			assign(paste("R_MovementTeff", "p1", position,combos2[indexpos], sep="_"),  ( 1-exp(-EBV/cEBV) ) * p_eff[paste(combos[indexpos])] * Teff )
			assign(paste("R_MovementTreg", "p1", position,combos2[indexpos], sep="_"),  ( 1-exp(-Teff/cTeff) ) * p_reg[paste(combos[indexpos])] * Treg )
			assign(paste("R_DACMovements","p1", position,combos2[indexpos], sep="_"),  ( 1-exp(-(Teff+Treg)/cTcell) ) * p_dac[paste(combos[indexpos])] * DAC )
			assign(paste("R_MovementEBV", "p1", position,combos2[indexpos], sep="_"),    p_ebv * EBV )
		}
		# end general transitions

	}

	####

	##Begin parameter rates

	TeffKillsODC1=TekODC
	TeffKillsODC2=TekODC
	TeffKillsODC3=TekODC
	TeffKillsODC4=TekODC

	TeffDup_Asym= Te2
	TeffDup_Sym=Te2

	TeffActivation= TeE
	MemActivation=MemE
	NKkillsTreg=NKkillT
	NKkillsTeff=NKkillT

	#Places
	EBV_P1=y[1]
	EBV_P2=y[2]
	EBV_P3=y[3]
	Teff_P1=y[4]
	Teff_P2=y[5]
	Teff_P3=y[6]
	Treg_P1=y[7]
	Treg_P2=y[8]
	Treg_P3=y[9]
	ODC_L1_P1=y[10]
	ODC_L1_P2=y[11]
	ODC_L1_P3=y[12]
	ODC_L2_P1=y[13]
	ODC_L2_P2=y[14]
	ODC_L2_P3=y[15]
	ODC_L3_P1=y[16]
	ODC_L3_P2=y[17]
	ODC_L3_P3=y[18]
	ODC_L4_P1=y[19]
	ODC_L4_P2=y[20]
	ODC_L4_P3=y[21]
	ODC_L5_P1=y[22]
	ODC_L5_P2=y[23]
	ODC_L5_P3=y[24]
	NK_P1=y[25]
	NK_P2=y[26]
	NK_P3=y[27]
	IL2_P1=y[28]
	IL2_P2=y[29]
	IL2_P3=y[30]
	DAC_P1=y[31]
	DAC_P2=y[32]
	DAC_P3=y[33]
	Resting_Teff_P1=y[34]
	Resting_Teff_P2=y[35]
	Resting_Teff_P3=y[36]
	Resting_Treg_P1=y[37]
	Resting_Treg_P2=y[38]
	Resting_Treg_P3=y[39]
	EffectorMemory=y[40]
	Resting_Treg_temp_P1=y[41]
	Resting_Treg_temp_P2=y[42]
	Resting_Treg_temp_P3=y[43]
	Resting_Teff_temp_P1=y[44]
	Resting_Teff_temp_P2=y[45]
	Resting_Teff_temp_P3=y[46]
	NK_temp_P1=y[47]
	NK_temp_P2=y[48]
	NK_temp_P3=y[49]



	#####  MASS ACTION:
	R_FromTimoEFF_p1_P1=prod(c( Resting_Teff_temp_P1^1))
	R_FromTimoEFF_p1_P2=prod(c( Resting_Teff_temp_P2^1))
	R_FromTimoEFF_p1_P3=prod(c( Resting_Teff_temp_P3^1))
	R_FromTimoREG_p1_P1=prod(c( Resting_Treg_temp_P1^1))
	R_FromTimoREG_p1_P2=prod(c( Resting_Treg_temp_P2^1))
	R_FromTimoREG_p1_P3=prod(c( Resting_Treg_temp_P3^1))
	R_Remyelinization_l1_L2_p1_P1=prod(c( ODC_L2_P1^1))
	R_Remyelinization_l1_L2_p1_P2=prod(c( ODC_L2_P2^1))
	R_Remyelinization_l1_L2_p1_P3=prod(c( ODC_L2_P3^1))
	R_Remyelinization_l1_L3_p1_P1=prod(c( ODC_L3_P1^1))
	R_Remyelinization_l1_L3_p1_P2=prod(c( ODC_L3_P2^1))
	R_Remyelinization_l1_L3_p1_P3=prod(c( ODC_L3_P3^1))
	R_Remyelinization_l1_L4_p1_P1=prod(c( ODC_L4_P1^1))
	R_Remyelinization_l1_L4_p1_P2=prod(c( ODC_L4_P2^1))
	R_Remyelinization_l1_L4_p1_P3=prod(c( ODC_L4_P3^1))
	R_TregDeath_p1_P1=prod(c( Treg_P1^1))
	R_TregDeath_p1_P2=prod(c( Treg_P2^1))
	R_TregDeath_p1_P3=prod(c( Treg_P3^1))
	R_TeffDeath_p1_P1=prod(c( Teff_P1^1))
	R_TeffDeath_p1_P2=prod(c( Teff_P2^1))
	R_TeffDeath_p1_P3=prod(c( Teff_P3^1))
	R_DACDegradation_p1_P1=prod(c( DAC_P1^1))
	R_DACDegradation_p1_P2=prod(c( DAC_P2^1))
	R_DACDegradation_p1_P3=prod(c( DAC_P3^1))
	R_NKDegradation_p1_P1=prod(c( NK_P1^1))
	R_NKDegradation_p1_P2=prod(c( NK_P2^1))
	R_NKDegradation_p1_P3=prod(c( NK_P3^1))
	R_NKarrive_p1_P1=prod(c( NK_temp_P1^1))
	R_NKarrive_p1_P2=prod(c( NK_temp_P2^1))
	R_NKarrive_p1_P3=prod(c( NK_temp_P3^1))

	## ODEs
	dEBV_P1 = + MovementEBV * ( 1*( 23*R_MovementEBV_p1_P3_p2_P1 ) + 1*( 1*R_MovementEBV_p1_P1_p2_P1 ) + 1*( 2*R_MovementEBV_p1_P2_p2_P1 ) )  - TeffkillsEBV * ( 1*( 1*R_TeffkillsEBV_p1_P1 ) )  - MovementEBV * ( 1*( 2*R_MovementEBV_p1_P1_p2_P2 ) + 1*( 1*R_MovementEBV_p1_P1_p2_P1 ) + 1*( 23*R_MovementEBV_p1_P1_p2_P3 ) )


	dEBV_P2 = + MovementEBV * ( 1*( 2*R_MovementEBV_p1_P1_p2_P2 ) + 1*( 1*R_MovementEBV_p1_P2_p2_P2 ) + 1*( 23*R_MovementEBV_p1_P3_p2_P2 ) )  - TeffkillsEBV * ( 1*( 1*R_TeffkillsEBV_p1_P2 ) )  - MovementEBV * ( 1*( 2*R_MovementEBV_p1_P2_p2_P1 ) + 1*( 23*R_MovementEBV_p1_P2_p2_P3 ) + 1*( 1*R_MovementEBV_p1_P2_p2_P2 ) )


	dEBV_P3 = + MovementEBV * ( 1*( 2*R_MovementEBV_p1_P1_p2_P3 ) + 1*( 22*R_MovementEBV_p1_P3_p2_P3 ) + 1*( 2*R_MovementEBV_p1_P2_p2_P3 ) )  - TeffkillsEBV * ( 1*( 1*R_TeffkillsEBV_p1_P3 ) )  - MovementEBV * ( 1*( 2*R_MovementEBV_p1_P3_p2_P2 ) + 1*( 22*R_MovementEBV_p1_P3_p2_P3 ) + 1*( 2*R_MovementEBV_p1_P3_p2_P1 ) )


	dTeff_P1 = + TeffActivation * ( 1*( 1*R_TeffActivation_p1_P1 ) )  + MemActivation * ( 1*( 1*R_MemActivation_p1_P1 ) )  + TeffDup_Sym * ( 1*( 1*R_TeffDup_Sym_p1_P1 ) )  + MovementTeff * ( 1*( 23*R_MovementTeff_p1_P3_p2_P1 ) + 1*( 1*R_MovementTeff_p1_P1_p2_P1 ) + 1*( 2*R_MovementTeff_p1_P2_p2_P1 ) )  - TeffDeath * ( 1*( 1*R_TeffDeath_p1_P1 ) )  - TeffkillsEBV * ( 1*( 1*R_TeffkillsEBV_p1_P1 ) )  - TregKillsTeff * ( 1*( 1*R_TregKillsTeff_p1_P1 ) )  - NKkillsTeff * ( 1*( 1*R_NKkillsTeff_p1_P1 ) )  - MovementTeff * ( 1*( 2*R_MovementTeff_p1_P1_p2_P2 ) + 1*( 1*R_MovementTeff_p1_P1_p2_P1 ) + 1*( 23*R_MovementTeff_p1_P1_p2_P3 ) )

	#cat("\ntime:",t,"\ndTeff_P1:",dTeff_P1,"\nTeffDup_Asym: ", TeffDup_Asym ," and ", ( 1*( R_TeffDup_Asym_p1_P1 ) ),"\nMemActivation:", MemActivation * ( 1*( R_MemActivation_p1_P1 ) ),"\nTeffDup_Sym", TeffDup_Sym," and ", ( 1*( 1*R_TeffDup_Sym_p1_P1 ) ),"\nTregKillsTeff",TregKillsTeff * ( 1*( 1*R_TregKillsTeff_p1_P1 ) ),"\nTeffkillsEBV",TeffkillsEBV * ( 1*( 1*R_TeffkillsEBV_p1_P1 ) ),"\n")
	
	
	dTeff_P2 = + TeffActivation * ( 1*( 1*R_TeffActivation_p1_P2 ) )  + MemActivation * ( 1*( 1*R_MemActivation_p1_P2 ) )  + TeffDup_Sym * ( 1*( 1*R_TeffDup_Sym_p1_P2 ) )  + MovementTeff * ( 1*( 2*R_MovementTeff_p1_P1_p2_P2 ) + 1*( 1*R_MovementTeff_p1_P2_p2_P2 ) + 1*( 23*R_MovementTeff_p1_P3_p2_P2 ) )  - TeffDeath * ( 1*( 1*R_TeffDeath_p1_P2 ) )  - TeffkillsEBV * ( 1*( 1*R_TeffkillsEBV_p1_P2 ) )  - TregKillsTeff * ( 1*( 1*R_TregKillsTeff_p1_P2 ) )  - NKkillsTeff * ( 1*( 1*R_NKkillsTeff_p1_P2 ) )  - MovementTeff * ( 1*( 2*R_MovementTeff_p1_P2_p2_P1 ) + 1*( 23*R_MovementTeff_p1_P2_p2_P3 ) + 1*( 1*R_MovementTeff_p1_P2_p2_P2 ) )

	#dTeff_P3 = + TeffActivation * ( 1*( 1*R_TeffActivation_p1_P3 ) )  + MemActivation * ( 1*( 1*R_MemActivation_p1_P3 ) )  + TeffDup_Sym * ( 1*( 1*R_TeffDup_Sym_p1_P3 ) + 2*( 1*R_TeffDup_Sym_p1_P2 ) + 2*( 1*R_TeffDup_Sym_p1_P1 ) )  + MovementTeff * ( 1*( 2*R_MovementTeff_p1_P1_p2_P3 ) + 1*( 22*R_MovementTeff_p1_P3_p2_P3 ) + 1*( 2*R_MovementTeff_p1_P2_p2_P3 ) )  - TeffDeath * ( 1*( 1*R_TeffDeath_p1_P3 ) )  - TeffkillsEBV * ( 1*( 1*R_TeffkillsEBV_p1_P3 ) )  - TregKillsTeff * ( 1*( 1*R_TregKillsTeff_p1_P3 ) )  - NKkillsTeff * ( 1*( 1*R_NKkillsTeff_p1_P3 ) )  - MovementTeff * ( 1*( 2*R_MovementTeff_p1_P3_p2_P2 ) + 1*( 22*R_MovementTeff_p1_P3_p2_P3 ) + 1*( 2*R_MovementTeff_p1_P3_p2_P1 ) )
	
	dTeff_P3 = + TeffActivation * ( 1*( 1*R_TeffActivation_p1_P3 ) )  + MemActivation * ( 1*( 1*R_MemActivation_p1_P3 ) )  + TeffDup_Sym * 1*( 1*R_TeffDup_Sym_p1_P3 ) + MovementTeff * ( 1*( 2*R_MovementTeff_p1_P1_p2_P3 ) + 1*( 22*R_MovementTeff_p1_P3_p2_P3 ) + 1*( 2*R_MovementTeff_p1_P2_p2_P3 ) )  - TeffDeath * ( 1*( 1*R_TeffDeath_p1_P3 ) )  - TeffkillsEBV * ( 1*( 1*R_TeffkillsEBV_p1_P3 ) )  - TregKillsTeff * ( 1*( 1*R_TregKillsTeff_p1_P3 ) )  - NKkillsTeff * ( 1*( 1*R_NKkillsTeff_p1_P3 ) )  - MovementTeff * ( 1*( 2*R_MovementTeff_p1_P3_p2_P2 ) + 1*( 22*R_MovementTeff_p1_P3_p2_P3 ) + 1*( 2*R_MovementTeff_p1_P3_p2_P1 ) )
	
	dTreg_P1 = + TregActivation * ( 1*( 1*R_TregActivation_p1_P1 ) )  + TregDup * ( 1*( 1*R_TregDup_p1_P1 ) )  + MovementTreg * ( 1*( 23*R_MovementTreg_p1_P3_p2_P1 ) + 1*( 1*R_MovementTreg_p1_P1_p2_P1 ) + 1*( 2*R_MovementTreg_p1_P2_p2_P1 ) )  - TregDeath * ( 1*( 1*R_TregDeath_p1_P1 ) )  - NKkillsTreg * ( 1*( 1*R_NKkillsTreg_p1_P1 ) )  - MovementTreg * ( 1*( 2*R_MovementTreg_p1_P1_p2_P2 ) + 1*( 1*R_MovementTreg_p1_P1_p2_P1 ) + 1*( 23*R_MovementTreg_p1_P1_p2_P3 ) )


	dTreg_P2 = + TregActivation * ( 1*( 1*R_TregActivation_p1_P2 ) )  + TregDup * ( 1*( 1*R_TregDup_p1_P2 ) )  + MovementTreg * ( 1*( 2*R_MovementTreg_p1_P1_p2_P2 ) + 1*( 1*R_MovementTreg_p1_P2_p2_P2 ) + 1*( 23*R_MovementTreg_p1_P3_p2_P2 ) )  - TregDeath * ( 1*( 1*R_TregDeath_p1_P2 ) )  - NKkillsTreg * ( 1*( 1*R_NKkillsTreg_p1_P2 ) )  - MovementTreg * ( 1*( 2*R_MovementTreg_p1_P2_p2_P1 ) + 1*( 23*R_MovementTreg_p1_P2_p2_P3 ) + 1*( 1*R_MovementTreg_p1_P2_p2_P2 ) )


	dTreg_P3 = + TregActivation * ( 1*( 1*R_TregActivation_p1_P3 ) )  + TregDup * ( 1*( 1*R_TregDup_p1_P3 ) )  + MovementTreg * ( 1*( 2*R_MovementTreg_p1_P1_p2_P3 ) + 1*( 22*R_MovementTreg_p1_P3_p2_P3 ) + 1*( 2*R_MovementTreg_p1_P2_p2_P3 ) )  - TregDeath * ( 1*( 1*R_TregDeath_p1_P3 ) )  - NKkillsTreg * ( 1*( 1*R_NKkillsTreg_p1_P3 ) )  - MovementTreg * ( 1*( 2*R_MovementTreg_p1_P3_p2_P2 ) + 1*( 22*R_MovementTreg_p1_P3_p2_P3 ) + 1*( 2*R_MovementTreg_p1_P3_p2_P1 ) )


	dODC_L1_P1 = + TeffKillsODC1 * ( 1*( 1*R_TeffKillsODC1_l1_L1_p1_P1 ) )

	dODC_L1_P2 = + TeffKillsODC1 * ( 1*( 1*R_TeffKillsODC1_l1_L1_p1_P2 ) )

	dODC_L1_P3 = + TeffKillsODC1 * ( 1*( 1*R_TeffKillsODC1_l1_L1_p1_P3 ) )

	dODC_L2_P1 = + TeffKillsODC2 * ( 1*( 1*R_TeffKillsODC2_l1_L2_p1_P1 ) )  - TeffKillsODC1 * ( 1*( 1*R_TeffKillsODC1_l1_L1_p1_P1 ) )  - Remyelinization * ( 1*( 1*R_Remyelinization_l1_L2_p1_P1 ) )


	dODC_L2_P2 = + TeffKillsODC2 * ( 1*( 1*R_TeffKillsODC2_l1_L2_p1_P2 ) )  - TeffKillsODC1 * ( 1*( 1*R_TeffKillsODC1_l1_L1_p1_P2 ) )  - Remyelinization * ( 1*( 1*R_Remyelinization_l1_L2_p1_P2 ) )


	dODC_L2_P3 = + TeffKillsODC2 * ( 1*( 1*R_TeffKillsODC2_l1_L2_p1_P3 ) )  - TeffKillsODC1 * ( 1*( 1*R_TeffKillsODC1_l1_L1_p1_P3 ) )  - Remyelinization * ( 1*( 1*R_Remyelinization_l1_L2_p1_P3 ) )


	dODC_L3_P1 = + TeffKillsODC3 * ( 1*( 1*R_TeffKillsODC3_l1_L3_p1_P1 ) )  + Remyelinization * ( 1*( 1*R_Remyelinization_l1_L2_p1_P1 ) )  - TeffKillsODC2 * ( 1*( 1*R_TeffKillsODC2_l1_L2_p1_P1 ) )  - Remyelinization * ( 1*( 1*R_Remyelinization_l1_L3_p1_P1 ) )


	dODC_L3_P2 = + TeffKillsODC3 * ( 1*( 1*R_TeffKillsODC3_l1_L3_p1_P2 ) )  + Remyelinization * ( 1*( 1*R_Remyelinization_l1_L2_p1_P2 ) )  - TeffKillsODC2 * ( 1*( 1*R_TeffKillsODC2_l1_L2_p1_P2 ) )  - Remyelinization * ( 1*( 1*R_Remyelinization_l1_L3_p1_P2 ) )


	dODC_L3_P3 = + TeffKillsODC3 * ( 1*( 1*R_TeffKillsODC3_l1_L3_p1_P3 ) )  + Remyelinization * ( 1*( 1*R_Remyelinization_l1_L2_p1_P3 ) )  - TeffKillsODC2 * ( 1*( 1*R_TeffKillsODC2_l1_L2_p1_P3 ) )  - Remyelinization * ( 1*( 1*R_Remyelinization_l1_L3_p1_P3 ) )


	dODC_L4_P1 = + TeffKillsODC4 * ( 1*( 1*R_TeffKillsODC4_l1_L4_p1_P1 ) )  + Remyelinization * ( 1*( 1*R_Remyelinization_l1_L3_p1_P1 ) )  - TeffKillsODC3 * ( 1*( 1*R_TeffKillsODC3_l1_L3_p1_P1 ) )  - Remyelinization * ( 1*( 1*R_Remyelinization_l1_L4_p1_P1 ) )


	dODC_L4_P2 = + TeffKillsODC4 * ( 1*( 1*R_TeffKillsODC4_l1_L4_p1_P2 ) )  + Remyelinization * ( 1*( 1*R_Remyelinization_l1_L3_p1_P2 ) )  - TeffKillsODC3 * ( 1*( 1*R_TeffKillsODC3_l1_L3_p1_P2 ) )  - Remyelinization * ( 1*( 1*R_Remyelinization_l1_L4_p1_P2 ) )


	dODC_L4_P3 = + TeffKillsODC4 * ( 1*( 1*R_TeffKillsODC4_l1_L4_p1_P3 ) )  + Remyelinization * ( 1*( 1*R_Remyelinization_l1_L3_p1_P3 ) )  - TeffKillsODC3 * ( 1*( 1*R_TeffKillsODC3_l1_L3_p1_P3 ) )  - Remyelinization * ( 1*( 1*R_Remyelinization_l1_L4_p1_P3 ) )


	dODC_L5_P1 = + Remyelinization * ( 1*( 1*R_Remyelinization_l1_L4_p1_P1 ) )  - TeffKillsODC4 * ( 1*( 1*R_TeffKillsODC4_l1_L4_p1_P1 ) )


	dODC_L5_P2 = + Remyelinization * ( 1*( 1*R_Remyelinization_l1_L4_p1_P2 ) )  - TeffKillsODC4 * ( 1*( 1*R_TeffKillsODC4_l1_L4_p1_P2 ) )


	dODC_L5_P3 = + Remyelinization * ( 1*( 1*R_Remyelinization_l1_L4_p1_P3 ) )  - TeffKillsODC4 * ( 1*( 1*R_TeffKillsODC4_l1_L4_p1_P3 ) )


	dNK_P1 = + NKdup * ( 1*( 1*R_NKdup_p1_P1 ) )  + NKarrive * ( 1*( 1*R_NKarrive_p1_P1 ) )  - NKkillsTreg * ( 1*( 1*R_NKkillsTreg_p1_P1 ) )  - NKkillsTeff * ( 1*( 1*R_NKkillsTeff_p1_P1 ) )  - NKDegradation * ( 1*( 1*R_NKDegradation_p1_P1 ) )

	#cat("\ntime:",t,"\ndNK_P1:",dNK_P1,"\nNKdup: ",NKdup* ( 1*( 1*R_NKdup_p1_P1 ) ),"\nNKarrive:", NKarrive * ( 1*( 1*R_NKarrive_p1_P1 ) ),"\nNKkillsTreg: ", NKkillsTreg * ( 1*( 1*R_NKkillsTreg_p1_P1 ) ),"\nNKkillsTeff", NKkillsTeff * ( 1*( 1*R_NKkillsTeff_p1_P1 ) ),"\nNKDegradation",NKDegradation * ( 1*( 1*R_NKDegradation_p1_P1 ) ) ,"\n")
	

	dNK_P2 = + NKdup * ( 1*( 1*R_NKdup_p1_P2 ) )  + NKarrive * ( 1*( 1*R_NKarrive_p1_P2 ) )  - NKkillsTreg * ( 1*( 1*R_NKkillsTreg_p1_P2 ) )  - NKkillsTeff * ( 1*( 1*R_NKkillsTeff_p1_P2 ) )  - NKDegradation * ( 1*( 1*R_NKDegradation_p1_P2 ) )


	dNK_P3 = + NKdup * ( 1*( 1*R_NKdup_p1_P3 ) )  + NKarrive * ( 1*( 1*R_NKarrive_p1_P3 ) )  - NKkillsTreg * ( 1*( 1*R_NKkillsTreg_p1_P3 ) )  - NKkillsTeff * ( 1*( 1*R_NKkillsTeff_p1_P3 ) )  - NKDegradation * ( 1*( 1*R_NKDegradation_p1_P3 ) )


	dIL2_P1 = + TeffActivation * ( 1*( 1*R_TeffActivation_p1_P1 ) )  + MemActivation * ( 1*( 1*R_MemActivation_p1_P1 ) )  - TeffDup_Asym * ( 1*( 1*R_TeffDup_Asym_p1_P1 ) )  - TregDup * ( 1*( 1*R_TregDup_p1_P1 ) )  - NKdup * ( 1*( 1*R_NKdup_p1_P1 ) )  - TeffDup_Sym * ( 1*( 1*R_TeffDup_Sym_p1_P1 ) )


	dIL2_P2 = + TeffActivation * ( 1*( 1*R_TeffActivation_p1_P2 ) )  + MemActivation * ( 1*( 1*R_MemActivation_p1_P2 ) )  - TeffDup_Asym * ( 1*( 1*R_TeffDup_Asym_p1_P2 ) )  - TregDup * ( 1*( 1*R_TregDup_p1_P2 ) )  - NKdup * ( 1*( 1*R_NKdup_p1_P2 ) )  - TeffDup_Sym * ( 1*( 1*R_TeffDup_Sym_p1_P2 ) )


	dIL2_P3 = + TeffActivation * ( 1*( 1*R_TeffActivation_p1_P3 ) )  + MemActivation * ( 1*( 1*R_MemActivation_p1_P3 ) )  - TeffDup_Asym * ( 1*( 1*R_TeffDup_Asym_p1_P3 ) )  - TregDup * ( 1*( 1*R_TregDup_p1_P3 ) )  - NKdup * ( 1*( 1*R_NKdup_p1_P3 ) )  - TeffDup_Sym * ( 1*( 1*R_TeffDup_Sym_p1_P3 ) )


	dDAC_P1 = + DACMovements * ( 1*( 23*R_DACMovements_p1_P3_p2_P1 ) + 1*( 1*R_DACMovements_p1_P1_p2_P1 ) + 1*( 2*R_DACMovements_p1_P2_p2_P1 ) )  - DACDegradation * ( 1*( 1*R_DACDegradation_p1_P1 ) )  - DACMovements * ( 1*( 2*R_DACMovements_p1_P1_p2_P2 ) + 1*( 1*R_DACMovements_p1_P1_p2_P1 ) + 1*( 23*R_DACMovements_p1_P1_p2_P3 ) )


	dDAC_P2 = + DACMovements * ( 1*( 2*R_DACMovements_p1_P1_p2_P2 ) + 1*( 1*R_DACMovements_p1_P2_p2_P2 ) + 1*( 23*R_DACMovements_p1_P3_p2_P2 ) )  - DACDegradation * ( 1*( 1*R_DACDegradation_p1_P2 ) )  - DACMovements * ( 1*( 2*R_DACMovements_p1_P2_p2_P1 ) + 1*( 23*R_DACMovements_p1_P2_p2_P3 ) + 1*( 1*R_DACMovements_p1_P2_p2_P2 ) )


	dDAC_P3 = + DACMovements * ( 1*( 2*R_DACMovements_p1_P1_p2_P3 ) + 1*( 22*R_DACMovements_p1_P3_p2_P3 ) + 1*( 2*R_DACMovements_p1_P2_p2_P3 ) )  - DACDegradation * ( 1*( 1*R_DACDegradation_p1_P3 ) )  - DACMovements * ( 1*( 2*R_DACMovements_p1_P3_p2_P2 ) + 1*( 22*R_DACMovements_p1_P3_p2_P3 ) + 1*( 2*R_DACMovements_p1_P3_p2_P1 ) )


	dResting_Teff_P1 = + FromTimoEFF * ( 1*( 1*R_FromTimoEFF_p1_P1 ) )  - TeffActivation * ( 1*( 1*R_TeffActivation_p1_P1 ) )


	dResting_Teff_P2 = + FromTimoEFF * ( 1*( 1*R_FromTimoEFF_p1_P2 ) )  - TeffActivation * ( 1*( 1*R_TeffActivation_p1_P2 ) )


	dResting_Teff_P3 = + FromTimoEFF * ( 1*( 1*R_FromTimoEFF_p1_P3 ) )  - TeffActivation * ( 1*( 1*R_TeffActivation_p1_P3 ) )


	dResting_Treg_P1 = + FromTimoREG * ( 1*( 1*R_FromTimoREG_p1_P1 ) )  - TregActivation * ( 1*( 1*R_TregActivation_p1_P1 ) )


	dResting_Treg_P2 = + FromTimoREG * ( 1*( 1*R_FromTimoREG_p1_P2 ) )  - TregActivation * ( 1*( 1*R_TregActivation_p1_P2 ) )


	dResting_Treg_P3 = + FromTimoREG * ( 1*( 1*R_FromTimoREG_p1_P3 ) )  - TregActivation * ( 1*( 1*R_TregActivation_p1_P3 ) )


	dEffectorMemory = + TeffDup_Asym * ( 1*( 2*R_TeffDup_Asym_p1_P1 ) )  + TeffDup_Asym * ( 1*( 2*R_TeffDup_Asym_p1_P2 ) )+ TeffDup_Asym * ( 1*( 23*R_TeffDup_Asym_p1_P3 ) ) - MemActivation * ( 1*( 2*R_MemActivation_p1_P1 ) ) - MemActivation * ( 1*( 2*R_MemActivation_p1_P2 ) )- MemActivation * ( 1*( 23*R_MemActivation_p1_P3 ) )


	dResting_Treg_temp_P1 = + TregActivation * ( 1*( 1*R_TregActivation_p1_P1 ) )  - FromTimoREG * ( 1*( 1*R_FromTimoREG_p1_P1 ) )


	dResting_Treg_temp_P2 = + TregActivation * ( 1*( 1*R_TregActivation_p1_P2 ) )  - FromTimoREG * ( 1*( 1*R_FromTimoREG_p1_P2 ) )


	dResting_Treg_temp_P3 = + TregActivation * ( 1*( 1*R_TregActivation_p1_P3 ) )  - FromTimoREG * ( 1*( 1*R_FromTimoREG_p1_P3 ) )


	dResting_Teff_temp_P1 = + TeffActivation * ( 1*( 1*R_TeffActivation_p1_P1 ) )  - FromTimoEFF * ( 1*( 1*R_FromTimoEFF_p1_P1 ) )


	dResting_Teff_temp_P2 = + TeffActivation * ( 1*( 1*R_TeffActivation_p1_P2 ) )  - FromTimoEFF * ( 1*( 1*R_FromTimoEFF_p1_P2 ) )


	dResting_Teff_temp_P3 = + TeffActivation * ( 1*( 1*R_TeffActivation_p1_P3 ) )  - FromTimoEFF * ( 1*( 1*R_FromTimoEFF_p1_P3 ) )


	dNK_temp_P1 = + NKkillsTreg * ( 1*( 1*R_NKkillsTreg_p1_P1 ) )  + NKkillsTeff * ( 1*( 1*R_NKkillsTeff_p1_P1 ) )  + NKDegradation * ( 1*( 1*R_NKDegradation_p1_P1 ) )  - NKarrive * ( 1*( 1*R_NKarrive_p1_P1 ) )-  NKdup * ( 1*( 1*R_NKdup_p1_P1 ) )


	dNK_temp_P2 = + NKkillsTreg * ( 1*( 1*R_NKkillsTreg_p1_P2 ) )  + NKkillsTeff * ( 1*( 1*R_NKkillsTeff_p1_P2 ) )  + NKDegradation * ( 1*( 1*R_NKDegradation_p1_P2 ) )  - NKarrive * ( 1*( 1*R_NKarrive_p1_P2 ) ) - NKdup * ( 1*( 1*R_NKdup_p1_P2 ) )


	dNK_temp_P3 = + NKkillsTreg * ( 1*( 1*R_NKkillsTreg_p1_P3 ) )  + NKkillsTeff * ( 1*( 1*R_NKkillsTeff_p1_P3 ) )  + NKDegradation * ( 1*( 1*R_NKDegradation_p1_P3 ) )  - NKarrive * ( 1*( 1*R_NKarrive_p1_P3 ) )-  NKdup * ( 1*( 1*R_NKdup_p1_P3 ) )



	list(c(dEBV_P1, dEBV_P2, dEBV_P3, dTeff_P1, dTeff_P2, dTeff_P3, dTreg_P1, dTreg_P2, dTreg_P3, dODC_L1_P1, dODC_L1_P2, dODC_L1_P3, dODC_L2_P1, dODC_L2_P2, dODC_L2_P3, dODC_L3_P1, dODC_L3_P2, dODC_L3_P3, dODC_L4_P1, dODC_L4_P2, dODC_L4_P3, dODC_L5_P1, dODC_L5_P2, dODC_L5_P3, dNK_P1, dNK_P2, dNK_P3, dIL2_P1, dIL2_P2, dIL2_P3, dDAC_P1, dDAC_P2, dDAC_P3, dResting_Teff_P1, dResting_Teff_P2, dResting_Teff_P3, dResting_Treg_P1, dResting_Treg_P2, dResting_Treg_P3, dEffectorMemory, dResting_Treg_temp_P1, dResting_Treg_temp_P2, dResting_Treg_temp_P3, dResting_Teff_temp_P1, dResting_Teff_temp_P2, dResting_Teff_temp_P3, dNK_temp_P1, dNK_temp_P2, dNK_temp_P3))
}

#Initial marking
yini<-c(y1=0, y2=0, y3=0, y4=0, y5=0, y6=0, y7=0, y8=0, y9=0, y10=0, y11=0, y12=0, y13=0, y14=0, y15=0, y16=0, y17=0, y18=0, y19=0, y20=0, y21=0, y22=500, y23=500, y24=500, y25=0, y26=0, y27=0, y28=0, y29=0, y30=0, y31=0, y32=0, y33=0, y34=0, y35=0, y36=0, y37=0, y38=0, y39=0, y40=0, y41=0, y42=0, y43=0, y44=0, y45=0, y46=0, y47=0, y48=0, y49=0)

y_names=c("EBV_P1", "EBV_P2", "EBV_P3", "Teff_P1", "Teff_P2", "Teff_P3", "Treg_P1", "Treg_P2", "Treg_P3", "ODC_L1_P1", "ODC_L1_P2", "ODC_L1_P3", "ODC_L2_P1", "ODC_L2_P2", "ODC_L2_P3", "ODC_L3_P1", "ODC_L3_P2", "ODC_L3_P3", "ODC_L4_P1", "ODC_L4_P2", "ODC_L4_P3", "ODC_L5_P1", "ODC_L5_P2", "ODC_L5_P3", "NK_P1", "NK_P2", "NK_P3", "IL2_P1", "IL2_P2", "IL2_P3", "DAC_P1", "DAC_P2", "DAC_P3", "Resting_Teff_P1", "Resting_Teff_P2", "Resting_Teff_P3", "Resting_Treg_P1", "Resting_Treg_P2", "Resting_Treg_P3", "EffectorMemory", "Resting_Treg_temp_P1", "Resting_Treg_temp_P2", "Resting_Treg_temp_P3", "Resting_Teff_temp_P1", "Resting_Teff_temp_P2", "Resting_Teff_temp_P3", "NK_temp_P1", "NK_temp_P2", "NK_temp_P3")

names(yini) = y_names

#### saving the positions of places
ebvPos<-c( which(y_names %in% grep("EBV_*", y_names, value=T)) )
teffPos<-c( which(y_names %in% grep("^Teff_*", y_names, value=T)) )
tregPos<-c( which(y_names %in% grep("^Treg_*", y_names, value=T)) )
odcPos<-c( which(y_names %in% grep("ODC_le1_*", y_names, value=T)))
il2Pos<-c( which(y_names %in% grep("IL2_*", y_names, value=T)))
tregPos<-c( which(y_names %in% grep("^Treg_*", y_names, value=T)) )
dacPos<-c( which(y_names %in% grep("DAC_*", y_names, value=T)) )
memPos<-c( which(y_names %in% grep("EffectorMemory", y_names, value=T)) )

yini[c("Resting_Treg_P1","Resting_Treg_P2","Resting_Treg_P3")] <- 63
yini[c("Resting_Teff_P1","Resting_Teff_P2","Resting_Teff_P3")] <- 1687
yini[c("NK_P1","NK_P2","NK_P3")] <- 375
yini[il2Pos] <- 1000

combos <-c("P1","P2","P3")

combos2<<- c("p2_P1","p2_P2","p2_P3")
