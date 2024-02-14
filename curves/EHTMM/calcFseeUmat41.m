function fsee = calcFseeUmat41(l_SEE,lSEE0,dUSEEnll,dUSEEl,dFSEE0)

l_SEE_nll = (1.0+dUSEEnll)*lSEE0;
v_SEE     = dUSEEnll/dUSEEl;
K_SEE_nl  = dFSEE0/((dUSEEnll*lSEE0)^v_SEE);
K_SEE_l   = dFSEE0/(dUSEEl*lSEE0);

%l_SEE = lp-lce;

if ( (l_SEE < l_SEE_nll) && (l_SEE > lSEE0) )
  fsee = K_SEE_nl*((l_SEE-lSEE0)^v_SEE);
elseif ((l_SEE > l_SEE_nll) && (l_SEE >lSEE0)) 
  fsee = dFSEE0+K_SEE_l*(l_SEE-l_SEE_nll);
else
  fsee = 0.0;
end
