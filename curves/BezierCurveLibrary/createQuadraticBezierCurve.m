function quadraticCurve = createQuadraticBezierCurve(higherOrderCurve)

assert(size(higherOrderCurve.xpts,1) > 3)

rootEPS = eps^0.5;

quadraticCurve = struct(...
    'xpts',zeros(3,size(higherOrderCurve.xpts,2)),...
    'ypts',zeros(3,size(higherOrderCurve.xpts,2)),...
    'xEnd',higherOrderCurve.xEnd,...
    'yEnd',higherOrderCurve.yEnd,...
    'dydxEnd',higherOrderCurve.dydxEnd,...
    'd2ydx2End',[0,0],...
    'integral',[]);

for i=1:1:size(higherOrderCurve.xpts,2)

    x0 = higherOrderCurve.xpts(1,i);
    x1 = higherOrderCurve.xpts(end,i);
    
    y0 = higherOrderCurve.ypts(1,i);
    y1 = higherOrderCurve.ypts(end,i);

    dydx0 = calcBezierYFcnXDerivative(x0,higherOrderCurve,1);

    dydx1 = calcBezierYFcnXDerivative(x1,higherOrderCurve,1);

    p0 = [x0,y0];
    p2 = [x1,y1];

    %1. Calculate the location where the two lines intersect
    % (x-x0)*dydx0 + y0 = (x-x1)*dydx1 + y1
    %   x*(dydx0-dydx1) = y1-y0-x1*dydx1+x0*dydx0
    %                 x = (y1-y0-x1*dydx1+x0*dydx0)/(dydx0-dydx1);
    
    xC = 0;
    yC1 = 0;
    yC2 = 0;
    
    if(abs(dydx0-dydx1) > rootEPS)
        xC = (y1-y0-x1*dydx1+x0*dydx0)/(dydx0-dydx1);    
        yC2 = (xC-x1)*dydx1 + y1;
        yC1=  yC2;
    else
        xC = (x1+x0)/2;
        yC1 = (xC-x0)*dydx0 + y0;
        yC2 = (xC-x1)*dydx1 + y1;
    end  
    if(abs(yC1-yC2) >= rootEPS)
        here=1;
    end
    assert(abs(yC1-yC2)<rootEPS);


    quadraticCurve.xpts(:,i) = [x0; xC;x1];
    quadraticCurve.ypts(:,i) = [y0;yC1;y1];

    if(i==1)
        quadraticCurve.dydxEnd(1,1)=dydx0;
    end
    if(i==size(higherOrderCurve.xpts,2))
        quadraticCurve.dydxEnd(1,2)=dydx1;
    end
   

end