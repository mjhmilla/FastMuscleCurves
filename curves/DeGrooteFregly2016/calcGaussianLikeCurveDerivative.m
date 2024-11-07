function y = calcGaussianLikeCurveDerivative(x,b1,b2,b3,b4)

t0 = (b2 - x);
t1 = (b3 + b4 * x);
t2 = t1*t1;
t3 = t2*t1;

y = (b1 * exp(-t0*t0 / (2 * t2 )) * t0 * (b3 + b2 * b4))/(t3);
