function fisom = calcFisomUmat41(lce,lceOpt,dWdes,nuCEdes, dWasc, nuCEasc)

fisom=0;
if (lce > lceOpt) 
    fisom = exp(-( abs( ((lce/lceOpt)-1.0)/dWdes ) )^nuCEdes);
else
    fisom = exp(-( abs( ((lce/lceOpt)-1.0)/dWasc) )^nuCEasc);
end

fisom = max(fisom,eps);
