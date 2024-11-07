function y = calcDeGrooteFregly2016ActiveForceLengthMultiplier(lceN,scale,der)

y = nan;

b11 = 0.8150671134243542;
b21 = 1.055033428970575;
b31 = 0.162384573599574;
b41 = 0.063303448465465;
b12 = 0.433004984392647;
b22 = 0.716775413397760;
b32 = -0.029947116970696;
b42 = 0.200356847296188;
b13 = 0.1;
b23 = 1.0;
b33 = 0.353553390593274; 
b43 = 0.0;

switch der
    case 0
        x = (lceN - 1.0) / scale + 1.0;
        y = calcGaussianLikeCurve(x, b11, b21, b31, b41) + ...
               calcGaussianLikeCurve(x, b12, b22, b32, b42) + ...
               calcGaussianLikeCurve(x, b13, b23, b33, b43);
    case 1

        x = (lceN - 1.0) / scale + 1.0;
        y = (1.0 / scale) * ...
               (calcGaussianLikeCurveDerivative(x, b11, b21, b31, b41) +...
                calcGaussianLikeCurveDerivative(x, b12, b22, b32, b42) +...
                calcGaussianLikeCurveDerivative(x, b13, b23, b33, b43));        
    otherwise
        assert(0, 'Error: der must be -1, 0, 1');
end
