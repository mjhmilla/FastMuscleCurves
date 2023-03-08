   
function fiberForceVelocityCurve =  createFiberForceVelocityCurve(...
                                        fvAtHalfVmax,...  
                                        fvAtLowEccVel,...                                        
                                        fvAtMaxEccVel,...
                                        eccCurviness,...
                                        flag_sharpEccentricTransition,...                                        
                                        flag_enableNumericallyNonZeroGradients,...
                                        smallNumericallyNonZeroValue,...
                                        smallNumericallyNonZeroSlope,...
                                        muscleName,...
                                        flag_usingOctave)

         


fiberForceVelocityCurve = [];
fiberForceVelocityCurve.name = sprintf('%s.%s',muscleName,...
                                       'fiberForceVelocityCurve');

vMax     = 1;                                   
vMaxC    = -vMax;   %since internally a concentric velocity is negative ...
                    %because the fiber is getting shorter.
vMaxE    =  vMax;
                  
assert( fvAtHalfVmax <= 0.45,...
        sprintf('%s: fvAtHalfVmax',...
        ' must be < 0.45',fiberForceVelocityCurve.name));
                                   

assert( fvAtMaxEccVel > 1.0, ...
  sprintf('%s: fvAtMaxEccVel must be greater than 1',fiberForceVelocityCurve.name));

assert( fvAtMaxEccVel > fvAtLowEccVel, ...
  sprintf('%s: fvAtMaxEccVel must be greater than fvAtLowEccVel',...
          fiberForceVelocityCurve.name));

dydxE = (fvAtMaxEccVel - fvAtLowEccVel)/vMaxE;

assert( (eccCurviness <= 1.0 && eccCurviness >= 0), ...    
  sprintf('%s: eccCurviness must be between 0 and 1',...
           fiberForceVelocityCurve.name));



%We are going to set these parameters so that they fit a Hill-type hyperbolic
%fv curve that has the same slope as dydxIso.
%
%First we compute the terms that are consistent with a Hill-type concentric 
%contraction. Starting from Hill's hyperbolic equation
%
% f(w)     = (fiso*b - a*w) / (b+w)
% df(w)/dw = [-(a)*(b+w) - (b*fiso-a*w)] / (b+w)^2
%
% since this is a normalized curve
%
% at w = vMaxC the numerator goes to 0
%
% (fiso*b - a*vMaxC) / (b + vMaxC) = 0
% (fiso*b - a*vMaxC) = 0;
%  b =  a*vMaxC/fiso;
%
% Subtituting this expression for b into the expression when
%  f(wHalf) = fvAtHalfVmax yields this expression for parameter a
%
%  a = fvAtHalfVmax*w*fiso ...
%      / (vMaxC*fvAtHalfVmax + fiso*vMaxC - fiso*w);
%
%%
fiso = 1;
w = 0.5*vMaxC;
a = -fvAtHalfVmax*w*fiso ...
    / (vMaxC*fvAtHalfVmax - fiso*vMaxC + fiso*w);
b =  a*vMaxC/fiso;

yCheck  = (b*fiso-a*w)/(b+w);
assert(abs(yCheck-fvAtHalfVmax) < sqrt(eps));

w = 0*vMaxC;
dydxIso =(-(a)*(b+w) - (b*fiso-a*w)) / ((b+w)*(b+w));

vNearMaxC = 0.8*vMaxC;

w         = vNearMaxC;
dydxNearC = (-(a)*(b+w) - (b*fiso-a*w)) / ((b+w)*(b+w));

yNearC = (b*fiso-a*w)/(b+w); 

vCExtrap = (-yNearC+dydxNearC*w)/dydxNearC;

assert( abs(dydxNearC) > abs(smallNumericallyNonZeroSlope) ...
         && abs(yNearC) > abs(smallNumericallyNonZeroValue) ...
         && vCExtrap > vMaxC, ...
  sprintf('%s: smallNumericallyNonZeroNumber must be less than %1.3e and %1.3e',...
  fiberForceVelocityCurve.name, dydxNearC, yCheck));


%Solve for the curviness parameter that results in the Bezier curve that 
%most closely matches Hill's concentric hyperbola on the concentric side. 

w      =  vNearMaxC;
xNearC =  w;
yNearC = (b*fiso-a*w)/(b+w);
xIso = 0;
yIso = 1.0;

cC = 0.5;
pts = calcQuinticBezierCornerControlPoints(xNearC,yNearC, dydxNearC, 0,...
                                             xIso,  yIso,   dydxIso, 0,...                                         
                                                                  cC);
curve.xpts = pts(:,1);
curve.ypts = pts(:,2);
curve.xEnd = [xNearC,xIso];
curve.yEnd = [yNearC,yIso];
curve.dydxEnd = [dydxNearC, dydxIso];

%Get the initial error between Hill and the Bezier curve for cC = 0.5.
nSample = 10;
f = 0;
for j=1:1:nSample
    w    = (j-1)*vMaxC/nSample;
    yHill = (b*fiso-a*w)/(b+w);        
    f = f + abs( calcBezierYFcnXDerivative(w, curve, 0) - yHill);
end

h = 0.25;
curveLeft.xpts = [];
curveLeft.ypts = [];
curveLeft.xEnd = [xNearC,xIso];
curveLeft.yEnd = [yNearC,yIso];
curveLeft.dydxEnd = [dydxNearC, dydxIso];

curveRight.xpts = [];
curveRight.ypts = [];
curveRight.xEnd = [xNearC,xIso];
curveRight.yEnd = [yNearC,yIso];
curveRight.dydxEnd = [dydxNearC, dydxIso];

%Use the bisection method to find the best curviness value
for i=1:1:10


  cCLeft = cC-h;
  ptsLeft =...
    calcQuinticBezierCornerControlPoints(xNearC,yNearC,dydxNearC, 0,...
                                         xIso,  yIso,  dydxIso, 0,...                                        
                                        cCLeft);
  curveLeft.xpts = ptsLeft(:,1);
  curveLeft.ypts = ptsLeft(:,2);
  
  cCRight = cC+h;
  ptsRight =...
    calcQuinticBezierCornerControlPoints(xNearC,yNearC,dydxNearC, 0,...
                                         xIso,  yIso,  dydxIso, 0,...                                       
                                       cCRight);
  curveRight.xpts = ptsRight(:,1);
  curveRight.ypts = ptsRight(:,2);

  %Compute the error at 10 points between -vMax and 0. 
  fLeft = 0;
  for j=1:1:nSample
    w     =  (j-1)*vMaxC/(nSample-1);
    yHill = (b*fiso-a*w)/(b+w);      
    fLeft = fLeft + abs( calcBezierYFcnXDerivative(w, curveLeft, 0) - yHill);
  end

  fRight = 0;
  for j=1:1:nSample
    w       = (j-1)*vMaxC/(nSample-1);
    yHill   = (b*fiso-a*w)/(b+w);     
    fRight  = fRight + abs( calcBezierYFcnXDerivative(w, curveRight, 0) - yHill);
  end
         
  %disp(sprintf('f: %e, fl: %e, fr: %e, cC: %e',f,fLeft,fRight,cC));  
  %Update the current solution
  if(abs(fLeft) < abs(f))
    f   = fLeft;
    cC  = cCLeft;
  end
  
  if(abs(fRight) < abs(f))
    f  = fRight;
    cC = cCRight; 
  end
  
  h     = h/2;

end


cE = scaleCurviness(eccCurviness);


xC    = vMaxC;
yC    = 0;
dydxC = 0;
                                        
if(flag_enableNumericallyNonZeroGradients==1)
  yC    = smallNumericallyNonZeroValue;
  dydxC = smallNumericallyNonZeroSlope;
end


concPts1 =...
    calcQuinticBezierCornerControlPoints(    xC,     yC,    dydxC,0,...
                                         xNearC, yNearC,dydxNearC,0, cC);
concPts2 =...
    calcQuinticBezierCornerControlPoints(xNearC,  yNearC,dydxNearC, 0, ...
                                         xIso,    yIso,  dydxIso, 0,  cC);


xE    = vMaxE;
yE    = fvAtMaxEccVel;
dydxE = (fvAtMaxEccVel-fvAtLowEccVel)/(vMaxE - 0);


if(flag_sharpEccentricTransition==1)

    xE0 = vMaxE*0.005;
    yE0 = yIso + 2*dydxIso*(xE0-xIso);
    dydxE0 = dydxIso*3;

    xE1    = vMaxE*0.1;
    yE1    = fvAtMaxEccVel+(xE1-xE)*dydxE;
    dydxE1 = dydxE;


    eccPts0  = ...
        calcQuinticBezierCornerControlPoints(xIso,     yIso, dydxIso, 0, ...
                                              xE0,      yE0,  dydxE0, 0, cE);
    eccPts1  = ...
        calcQuinticBezierCornerControlPoints(xE0,     yE0, dydxE0, 0, ...
                                             xE1,      yE1, dydxE1, 0, cE);
  
                                                                              
    xpts = [concPts1(:,1) concPts2(:,1)  eccPts0(:,1) eccPts1(:,1) ];
    ypts = [concPts1(:,2) concPts2(:,2)  eccPts0(:,2) eccPts1(:,2) ];
    
    fiberForceVelocityCurve.xpts = xpts;
    fiberForceVelocityCurve.ypts = ypts;

    xE=xE1;
    yE=yE1;
    dydxE=dydxE1;

else
    eccPts0  = ...
        calcQuinticBezierCornerControlPoints(xIso,     yIso, dydxIso, 0, ...
                                               xE,       yE,   dydxE, 0, cE);
                                                                              
    xpts = [concPts1(:,1) concPts2(:,1)  eccPts0(:,1)];
    ypts = [concPts1(:,2) concPts2(:,2)  eccPts0(:,2)];
    
    fiberForceVelocityCurve.xpts = xpts;
    fiberForceVelocityCurve.ypts = ypts;
end


fiberForceVelocityCurve.xEnd = [xC xE];
fiberForceVelocityCurve.yEnd = [yC yE];
fiberForceVelocityCurve.dydxEnd  = [dydxC dydxE];
fiberForceVelocityCurve.d2ydx2End= [0,0];
fiberForceVelocityCurve.integral = [];
