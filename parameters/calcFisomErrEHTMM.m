function errVal= calcFisomErrEHTMM(arg, fisomParams,targetPoints)

dWdes   = fisomParams.dWdes;
nuCEdes = fisomParams.nuCEdes;
dWasc   = fisomParams.dWasc;
nuCEasc = fisomParams.nuCEasc;
lceOpt  = fisomParams.lceOpt;

dWasc   = arg(1,1);
nuCEasc = arg(2,1);

lceN0Target = targetPoints.lceN0;
flN0Target  = targetPoints.flN0;
lceN1Target = targetPoints.lceN1;
flN1Target  = targetPoints.flN1;
lceN2Target = targetPoints.lceN2;
flN2Target  = targetPoints.flN2;


flN0  = calcFisomUmat41(lceN0Target*lceOpt,lceOpt,dWdes,nuCEdes, dWasc, nuCEasc);
err0  = flN0-flN0Target;

flN1  = calcFisomUmat41(lceN1Target*lceOpt,lceOpt,dWdes,nuCEdes, dWasc, nuCEasc);
err1  = flN1-flN1Target;

flN2  = calcFisomUmat41(lceN2Target*lceOpt,lceOpt,dWdes,nuCEdes, dWasc, nuCEasc);
err2  = flN2-flN2Target;

errVal = [err0;err1;err2];