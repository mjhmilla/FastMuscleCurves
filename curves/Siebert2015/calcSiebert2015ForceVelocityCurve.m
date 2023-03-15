function fvN = calcSiebert2015ForceVelocityCurve(velInLPS, fvParameters, normEccentricDamping)


fvN = 0;

if(velInLPS <= 0)
    if(abs(velInLPS) < fvParameters.vccmax)
        fvN = (fvParameters.vccmax-abs(velInLPS))...
             /(fvParameters.vccmax + abs(velInLPS)/fvParameters.curv);
    end

else
    %Siebert et al. does not describe what happens in the eccentric side
    fvN = 1. + normEccentricDamping*velInLPS/fvParameters.vccmax;    
end
