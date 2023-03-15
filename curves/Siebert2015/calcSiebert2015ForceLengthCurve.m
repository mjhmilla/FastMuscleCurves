function fpeN = calcSiebert2015ForceLengthCurve(lengthInMM, fpeParameters)

fpeN = 0;

dlpe = (lengthInMM-fpeParameters.lpec0);

if(dlpe > 0)
    fpeN = fpeParameters.k1*(exp(fpeParameters.k2*dlpe)-1);
    fpeN = fpeN/fpeParameters.Fim;
end


