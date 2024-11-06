function y = calcDeGrooteFregly2016ForceVelocityMultiplier(...
                    vceN,d1,d2,d3,d4, der)

y = nan;

switch der
    case 0
        tempV = d2 * vceN + d3;
        tempLogArg = tempV + sqrt(tempV*tempV + 1.0);
        y = d1 * log(tempLogArg) + d4;
    otherwise
        assert(0, 'Error: der must be -1, 0, 1');
end
