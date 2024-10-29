function errVal= calcFpeErrEHTMM(arg, fpeeParams,targetPoints)

lceOpt  = fpeeParams.lceOpt;
dWdes   = fpeeParams.dWdes;
nuCEdes = fpeeParams.nuCEdes;
Fmax    = fpeeParams.Fmax;
FPEE    = fpeeParams.FPEE;
LPEE0   = fpeeParams.LPEE0;
nuPEE   = fpeeParams.nuPEE;
lSEE0   = fpeeParams.lSEE0;
dUSEEnll= fpeeParams.dUSEEnll;
dUSEEl  = fpeeParams.dUSEEl;
dFSEE0  = fpeeParams.dFSEE0;

LPEE0=arg(1,1);
FPEE = arg(2,1);

lceN0Target=targetPoints.lceN0;
fpeN0Target=targetPoints.fpeN0;
fpeN1Target=targetPoints.fpeN1;
dfpeN1Target=targetPoints.dfpeN1;

%Subtract off the tendon strain
lt0         = calcFseeInverseUmat41(fpeN0Target,lSEE0,dUSEEnll,dUSEEl,dFSEE0);
dlt0        = lt0-lSEE0;
lce0Target  = lceN0Target*lceOpt - dlt0;

fpe0  = calcFpeeUmat41(lce0Target,lceOpt,dWdes,Fmax,FPEE,LPEE0,nuPEE); 
fpeN0 = fpe0/Fmax;
err0  = fpeN0-fpeN0Target;


lceFmax     = calcFpeeInverseUmat41(fpeN1Target*Fmax, lceOpt,dWdes,Fmax,FPEE,LPEE0,nuPEE);
dfpeFmax    = calcFpeeDerivativeUmat41(lceFmax, lceOpt,dWdes,Fmax,FPEE,LPEE0,nuPEE);

ltFmax      = calcFseeInverseUmat41(fpeN1Target*Fmax,lSEE0,dUSEEnll,dUSEEl,dFSEE0);
dftFmax     = calcFseeDerivativeUmat41(ltFmax,lSEE0,dUSEEnll,dUSEEl,dFSEE0); 

dfmtFmax = ((1/dfpeFmax)+(1/dftFmax))^(-1);

dfmtNFmax   = dfmtFmax*(lceOpt/Fmax);

err1    = dfmtNFmax-dfpeN1Target; 

errVal = [err0; err1];