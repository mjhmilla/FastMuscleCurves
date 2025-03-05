function tanhCoeffs = ...
    fitTanhSplineCoefficients(xk,yk,dydxk,yLimits,xAtIntYZero,xSample,ySample)

assert(size(xk,2)==1,'Error: xk, yk, dydxk must be n x 1');
assert(length(xk)==length(yk) && length(xk)==length(dydxk),...
    'Error: xk, yk, and dydxk must be the same length');
assert(size(yLimits,1)==1 && size(yLimits,2)==2,'Error:' )

n = length(xk)-1;


xS = ones(size(xk))*1;
dy = diff(yk).*0.05;

tanhCoeffs = zeros(n, 6);



for i=2:1:length(xk)

    j = i-1; %Segment counter


    x0 = xk(i-1,1);
    x1 = xk(i,1);

    y0 = yk(i-1,1);
    y1 = yk(i,1);

    dydx0 = dydxk(i-1,1);
    dydx1 = dydxk(i,1);

    %Calculate the location where the two lines intersect
    xC = 0;
    yC1 = 0;
    yC2 = 0;
    rootEPS = eps^0.5;
    if(abs(dydx0-dydx1) > rootEPS)
        xC = (y1-y0-x1*dydx1+x0*dydx0)/(dydx0-dydx1);    
        yC2 = (xC-x1)*dydx1 + y1;
        yC1=  yC2;
    else
        xC = (x1+x0)/2;
        yC1 = (xC-x0)*dydx0 + y0;
        yC2 = (xC-x1)*dydx1 + y1;
    end
    
    lambda = 0.;
    x0T = x0 + (xC-x0)*lambda;
    x1T = x1 + (xC-x1)*lambda;

    y0T = y0 + dydx0*(xC-x0T);
    y1T = y1 + dydx1*(xC-x1T);

    dydxT0 = 0;
    dydxT1 = 0;

    yNegInf = 0;
    yInf    = sign(dydx1)*inf;    
    
    xScale = xS(j,1);

    if(j==1)
        yNegInf = yLimits(1,1);
        dydxT0 = dydx0;
        dydxT1 = dydx1;
    else
        dydxT0 = 0;
        dydxT1 = dydx1-calcTanhSeriesDerivative(x1,tanhCoeffs(1:(j-1),:),1);
        
        yInf    = sign(dydxT1)*inf;        
    end

    xPoint = 0;
    yPoint = 0;
    if(j==1)
        xPoint = x1;
        yPoint = y1;
    else
        %Compensate for the error of the end point of the last curve
        xPoint = x1;
        yPoint = y1 - calcTanhSeriesDerivative(xPoint,tanhCoeffs(1:(j-1),:),0);
    end

    disp(j);
    fprintf('%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\n',...
        x0,x1,dydxT0,dydxT1,yNegInf,yInf,xScale,xPoint,yPoint);

    
    [A,B,C,D,E,F] = calcTanhSegmentCoefficientsUpd( ...
                        x0T,     ...
                        x1T,     ...
                        dydxT0,  ...
                        dydxT1,  ...
                        yNegInf,...
                        yInf,   ...
                        xScale, ...
                        xPoint, ...
                        yPoint, ...
                        xAtIntYZero)  ;  
    tanhCoeffs(j,:) = [A,B,C,D,E,F];

    x0P = x0;
    x1P = x1;
    y0P = y0;
    y1P = y1;
    dydx0P = dydx0;
    dydx1P = dydx1;

    flag_debug = 0;
    if(flag_debug==1)
        if(j==1)
            figDebug=figure;
        end
        xT = [xk(1,1):((x1-xk(1,1))/99):x1]';
        yT = zeros(size(xT));
        yS = zeros(size(xT));
        for z=1:1:length(xT)
            yT(z,1) = calcTanhSeriesDerivative(xT(z,1),tanhCoeffs(j,:),0);
            yS(z,1) = calcTanhSeriesDerivative(xT(z,1),tanhCoeffs(1:j,:),0);            
        end
        figure(figDebug);
        plot(xT,yS,'-k');
        hold on;
        plot(xT,yT,'-r');
        hold on;
        plot(xPoint,yPoint,'xk');
        hold on;
        plot(x0,y0,'ok');
        hold on;
        plot(x1,y1,'^k');        
        hold on;
        xlabel('X');
        ylabel('Y');
        box off;
        pause(0.1);
        clf(figDebug);

    end
end






