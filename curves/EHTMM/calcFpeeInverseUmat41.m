function lce = calcFpeeInverseUmat41(fpee, lceOpt,dWdes,Fmax,FPEE,LPEE0,nuPEE)


lce = LPEE0*lceOpt;

if(fpee > 0)
	KPEE = (FPEE*Fmax) / (lceOpt*(dWdes+1.0-LPEE0))^nuPEE;
	lce = (fpee/KPEE)^(1/nuPEE)+LPEE0*lceOpt;

	%fpee = KPEE*((lce-LPEE0)^nuPEE);
end