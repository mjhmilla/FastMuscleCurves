function ftN = calcSiebert2015TendonForceLengthCurve(lengthInMM, ftParameters)

ftN = 0;

if(lengthInMM > ftParameters.lsec0)
    dlsec   = (lengthInMM-ftParameters.lsec0)/ftParameters.lsec0;

    lsec1   = ftParameters.lsec0*ftParameters.lsec1_lsec0;
    dlsec1  = (lsec1)/ftParameters.lsec0;

    F1_Fim  = ftParameters.F1_Fim;

    if(dlsec < dlsec1)
        ksh = ftParameters.ksh;
        ftN = (F1_Fim/(exp(ksh)-1))*( exp(ksh*dlsec/dlsec1) -1 );
    else
%        k   = ftParameters.k;

        ksh = ftParameters.ksh;
        k = (F1_Fim/(exp(ksh)-1))*(ksh/dlsec1)*exp(ksh);
        ftN = F1_Fim + k*(dlsec-dlsec1);
    end

end


