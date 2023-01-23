function [A,B,C,D,E,F] = ...
    calcTanSegmentCoefficients(x0,x1,dydx0,dydx1,yNegInf,yInf, xAtIntYZero, xShift, xScale)

assert(isempty(yNegInf) || isempty(yInf), ...
    ['Error: this function has been formulated so that the function ',...
     'through x0 and x1 with slopes of dydx0 and dydx1. As a result, ',...
     'this function can be constrained to pass through yNegInf or yInf, but '...
     'not both. Set the yNegInf or yInf that you do not care about to zero. ']);

%B and C can be adjusted
B = x0 + (x1-x0)*0.5 + xShift;
C = ((x1-x0)*0.5)*(1/2)*xScale;

%A and D must take these values to meet the derivative constraints in the 
%limit


A = (dydx1-dydx0)*(1./pi);
D = dydx0 + (dydx1-dydx0)*0.5;

if(isempty(yNegInf) == 0)
    assert(abs(dydx0) < eps);
    %This means that yNegInf has a finite value, which means that dydx0=0
    %D = 0, and A = dydx1/pi.
    %
    % y = A*( argXB*atan(argXBC)-log(argXBC*argXBC+1)*0.5*C)+D*x + E
    %   = A*( argXB*atan(argXBC)-log(argXBC*argXBC+1)*0.5*C)+D*x + E
    %
    % as x->-inf term by term we have
    %   = A*( ...
    %        argXB*atan(argXBC) ...         : -> -pi*0.5*(x-B)/C
    %       -log(argXBC*argXBC+1)*0.5*C ... : -> -2*log( ((x-B)/C)*((x-B)/C) )
    %      )...         
    %      +D*x                             : 0, since D = 0
    %      + E                              :
    %
    % I'm short on time right now, so I'm just going to do this numerically
    x       = -1000.0*C + B;
    argXB   = x-B;
    argXBC  = argXB/C;    
    E = -(A*( argXB*atan(argXBC)-log(argXBC*argXBC+1)*0.5*C)+D*x);
end


if(isempty(yInf) == 0)
    assert(abs(dydx1) < eps);

    % I'm short on time right now, so I'm just going to do this numerically
    x       = 1000.0*C + B;
    argXB   = x-B;
    argXBC  = argXB/C;    
    E = -(A*( argXB*atan(argXBC)-log(argXBC*argXBC+1)*0.5*C)+D*x);    
end

% I'm short on time right now, so I'm just going to do this numerically
x = xAtIntYZero;

argXBC = (x-B)/C;

dC = 1/C;
C2 = C*C;
C3 = C2*C;
C4 = C3*C;
B2 = B*B;

x2 = x*x;

atanXBC = atan(argXBC);

inty = A*C*(...
          ((x2*0.5-B*x) * atanXBC ...
              -(...
                  (...
                    (-C4-B2*C2)*atanXBC...
                  )*(0.5*dC) ...
                  +(C2*x)*0.5 ...
                )*dC ...
          )*dC ...
         -( ...
            x*log(argXBC*argXBC+1) ...
            -(2*(...
                  (B*C2*log(x2-2*B*x+C2+B2))*0.5 ...
                  -C3*atanXBC+C2*x ...
                ) ...
              )*(dC*dC) ...
          )*0.5 ...
        )...
        +(D*x2)*0.5...
        +E*x;

F = -inty;








