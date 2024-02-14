function fpee = calcFpeeUmat41(lce, lceOpt,dWdes,Fmax,FPEE,LPEE0,nuPEE)


fpee=0;

if(lce > lceOpt*LPEE0)
	KPEE = (FPEE*Fmax) / (lceOpt*(dWdes+1.0-LPEE0))^nuPEE;
	fpee = KPEE*((lce-lceOpt*LPEE0)^nuPEE);
end