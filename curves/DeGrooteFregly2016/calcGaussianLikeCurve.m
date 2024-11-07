function y = calcGaussianLikeCurve(x,b1,b2,b3,b4)

t0 = (x - b2);
t1 = (b3 + b4 * x);
y = b1 * exp(-0.5 * t0*t0 / (t1*t1));