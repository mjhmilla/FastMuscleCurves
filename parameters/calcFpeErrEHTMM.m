function errVal= calcFpeErrEHTMM(arg, fpeeParams,targetPoints)

lceOpt  = fpeeParams.lceOpt;
dWdes   = fpeeParams.dWdes;
nuCEdes = fpeeParams.nuCEdes;
Fmax    = fpeeParams.Fmax;
FPEE    = fpeeParams.FPEE;
LPEE0   = fpeeParams.LPEE0;
nuPEE   = fpeeParams.nuPEE;

LPEE0=arg(1,1);
FPEE = arg(2,1);

lceN0Target=targetPoints.lceN0;
fpeN0Target=targetPoints.fpeN0;
lceN1Target=targetPoints.lceN1;
fpeN1Target=targetPoints.fpeN1;
dfpeN1Target=targetPoints.dfpeN1;

fpe0  = calcFpeeUmat41(lceN0Target*lceOpt, lceOpt,dWdes,Fmax,FPEE,LPEE0,nuPEE); 
fpeN0 = fpe0/Fmax;
err0  = fpeN0-fpeN0Target;

lceFmax     = calcFpeeInverseUmat41(fpeN1Target*Fmax, lceOpt,dWdes,Fmax,FPEE,LPEE0,nuPEE);
dfpeFmax    = calcFpeeDerivativeUmat41(lceFmax, lceOpt,dWdes,Fmax,FPEE,LPEE0,nuPEE);
dfpeNFmax   = dfpeFmax*(lceOpt/Fmax);

err1    = dfpeNFmax-dfpeN1Target; 

errVal = [err0; err1];