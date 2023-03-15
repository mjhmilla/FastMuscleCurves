function falN = calcSiebert2015ActiveForceLengthCurve(lengthInMM, falParameters)

falN    = nan;
lN      = lengthInMM/falParameters.lccopt;


if(lN < falParameters.l1N)
    falN = 0;

elseif(lN > falParameters.l1N  && lN <= falParameters.l2N)
    ld      = (lN-falParameters.l1N)/(falParameters.l2N-falParameters.l1N);
    falN    = ld*falParameters.fc;

elseif(lN > falParameters.l2N  && lN <= 1)
    ld      = (lN-falParameters.l2N)/(1-falParameters.l2N);
    falN    = falParameters.fc + ld*(1-falParameters.fc);

elseif(lN > 1                  && lN <= falParameters.l3N)
    falN    = 1;

elseif(lN > falParameters.l3N  && lN <= falParameters.l4N)
    ld      = (lN-falParameters.l3N)/(falParameters.l4N-falParameters.l3N);
    falN    = 1 - ld;

else
    falN = 0;    
end

