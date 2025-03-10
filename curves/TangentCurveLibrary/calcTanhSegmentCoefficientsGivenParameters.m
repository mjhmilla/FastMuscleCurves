function [A,B,C,D,E,F] = ...
    calcTanhSegmentCoefficientsGivenParameters(x0,x1,...
                                dydx0,dydx1,...
                                yNegInf,yInf, ...
                                xScale, xPoint,yPoint,...
                                xAtIntYZero,...
                                indexParameters, valueParameters)

%fprintf('%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\n',...
%    x0,x1,dydx0,dydx1,yNegInf,yInf,xScale,xPoint,yPoint);

%for i=1:1:length(indexParameters)
%    assert(indexParameters(i,1) == 1 || indexParameters(i,1) == 3,...
%        ['Error: Only A and C may be set directly.']);
%end

xShift=0;
xShiftWidth=(max(x1,xPoint)-min(x0,xPoint));

assert( (yPoint >= yNegInf && yPoint <= yInf && (yInf-yNegInf) > 0) ...
       || (yPoint <= yNegInf && yPoint >= yInf  && (yNegInf-yInf) > 0),...
       'Error: yPoint must be between yNegInf and yInf');

assert(isinf(yNegInf) || isinf(yInf), ...
    ['Error: this function has been formulated so that the function ',...
     'through x0 and x1 with slopes of dydx0 and dydx1. As a result, ',...
     'this function can be constrained to pass through yNegInf or yInf, but '...
     'not both. Set the yNegInf or yInf that you do not care about to zero. ']);

xShiftBest=0;
B = x0 + (x1-x0)*0.5+xShiftBest;
C = ((x1-x0)*0.5)*(1/2)*xScale;
A = -(dydx0-dydx1)*0.5;

for i=1:1:length(indexParameters)
    switch indexParameters(i,1)
        case 1
            A = valueParameters(i,1);
        case 2
            B = valueParameters(i,1);
        case 3
            C = valueParameters(i,1);            
    end
end

D = dydx1-A;
E=0;
if(~isinf(yNegInf))
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




if(~isinf(yInf))
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

yValue          = calcTanhSeriesDerivative(xPoint,[A,B,C,D,E,nan],0);
yErrBest        = abs(yValue-yPoint);
xShiftBest      = 0;
BBest = B;

for i=1:1:20
    %B and C can be adjusted
    %B = x0 + (x1-x0)*0.5;% + xShift;
        
    BLeft = BBest - xShiftWidth;
    BRight= BBest + xShiftWidth;

    yL          = calcTanhSeriesDerivative(xPoint,[A,BLeft,C,D,E,nan],0);
    yLErr       = abs(yL-yPoint);
    yR          = calcTanhSeriesDerivative(xPoint,[A,BRight,C,D,E,nan],0);
    yRErr       = abs(yR-yPoint);
    
    if(yLErr < yErrBest && yLErr <= yRErr)
        yErrBest=yLErr;
        BBest=BLeft;
    elseif(yRErr < yErrBest && yRErr < yLErr)
        yErrBest=yRErr;
        BBest=BRight;        
    end

    xShiftWidth=xShiftWidth*0.5;
end

B=BBest;

x=xAtIntYZero;
C2 = C*C;

inty= 0.5.*((A*C2).*polylog(2,-exp((2.*(B-x))./C)) ...
                     + A.*(B-x).*(B+C*log(4)-x) ...
                     + x.*(2.*(-B*D+C+E) + D.*x)...
                     );

F = -inty;








