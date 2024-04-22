function val = calcLinearlyExtrapolatedExponential(argX, coeff)

c=coeff.c;
k=coeff.k;
ls=coeff.ls;

val = (c*k).*log(exp( (argX-ls)./k ) +1);
