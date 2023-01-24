function [A,B,C,D,E,F] = ...
    calcTanhSegmentCoefficients(x0,x1,dydx0,dydx1,yNegInf,yInf, xShift, xScale, xAtIntYZero)

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
A = -(dydx0-dydx1)*0.5;
D = dydx1-A;

if(isempty(yNegInf) == 0)
    assert(abs(dydx0) < eps);
    %This means that yNegInf has a finite value, which means that dydx0=0
    %
    %y = D*(x-B) + A*C*( log( ( exp((x-B)/C)+exp(-(x-B)/C) )/2 ) + C + E
    %
    %In the limit, as x-> -inf, the term exp((x-B)/C) goes to zero and
    %we have
    %
    % y = D*(x-B) + A*C*( log( exp(-(x-B)/C  )/2      )  + C + E
    % y = D*(x-B) + A*C*( log( exp(-(x-B)/C) )-log(2) )  + C + E
    % y = D*(x-B) + A*C*(-((x-B)/C) -log(2)           )  + C + E
    %
    % Since we want y=yNegInf as x->-inf, we set E to
    % E = dydx0 - (D*(x-B) + A*C*(-((x-B)/C)-log(2)) + C )
    %   = dydx0 - (D*(x-B) + -A*C*((x-B)/C)- A*C*log(2) + C)
    %   = dydx0 - (D*(x-B) + -A*(x-B)- A*C*log(2) + C)
    % Since D = A this simplifies to
    %   = dydx0 - (- A*C*log(2) + C)
    E = dydx0 - ( -A*C*log(2) + C );
end


if(isempty(yInf) == 0)
    assert(abs(dydx1) < eps);
    %This means that yNegInf has a finite value, which means that dydx0=0
    
    %y = D*(x-B) + A*C*( log( ( exp((x-B)/C)+exp(-(x-B)/C) )/2 ) + C + E
    %
    %In the limit, as x-> inf, the term exp(-(x-B)/C) goes to zero and
    %we have
    %
    % y = D*(x-B) + A*C*( log( exp( (x-B)/C  )/2      )  + C + E
    % y = D*(x-B) + A*C*( log( exp( (x-B)/C) )-log(2) )  + C + E
    % y = D*(x-B) + A*C*( ((x-B)/C) -log(2)           )  + C + E
    %
    % Since we want y=yInf as x->inf, we set E to
    % E = dydx1 - (D*(x-B) + A*C*( (x-B)/C) -log(2)   ) + C )
    %   = dydx0 - (D*(x-B) + A*C*((x-B)/C)- A*C*log(2) + C)
    %   = dydx0 - (D*(x-B) + A*(x-B) - A*C*log(2) + C)
    % Since D = -A this simplifies to
    %   = dydx0 - (- A*C*log(2) + C)
    E = dydx1   - (-A*C*log(2) + C );
end

x=xAtIntYZero;
C2 = C*C;

inty= 0.5.*((A*C2).*polylog(2,-exp((2.*(B-x))./C)) ...
                     + A.*(B-x).*(B+C*log(4)-x) ...
                     + x.*(2.*(-B*D+C+E) + D.*x)...
                     );

F = -inty;








