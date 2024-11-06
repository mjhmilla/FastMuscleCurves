function y = calcDeGrooteFregly2016PassiveMultiplier(lceN,e0,kPE,lceNMin,der)

y = nan;

switch der
    case -1
        temp1 = exp(kPE * lceNMin / e0);
        denom = exp(kPE * (1.0 + 1.0 / e0)) - temp1;
        temp2 = kPE / e0 * lceN;
        y= (e0 / kPE * exp(temp2) - lceN * temp1) / denom;

    case 0
        offset  = exp(kPE * (lceNMin - 1.0) / e0);
        denom   = exp(kPE) - offset;
        y       = (exp(kPE*(lceN-1.0)/e0) - offset) / denom; 

    case 1

        offset = exp(kPE * (lceNMin - 1) / e0);

        y = (kPE*exp((kPE*(lceN-1))/e0))/(e0*(exp(kPE)-offset));

    otherwise
        assert(0, 'Error: der must be -1, 0, 1');
end
