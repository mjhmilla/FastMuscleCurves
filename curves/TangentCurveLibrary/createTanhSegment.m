function [tanhCoeffs,curveInputs] = createTanhSegment(tanhCoeffs, ...
                             indexOfSegment,...
                             indicesOfKnotPoints,...
                             xk, yk, dydxk, smoothness,...
                             yLimits,xAtIntYZero,...
                             curveInputs,...
                             flag_appendCoeffcients)

    i  = indicesOfKnotPoints(1,2);
    j = indexOfSegment; %Segment counter

    assert(diff(indicesOfKnotPoints) == 1,'The indices must differ by 1');
    assert((i-j) == 1, 'The segment must be 1 less than the first knot index');

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
    
    lambda = 0;
    x0T = x0 + (xC-x0)*lambda;
    x1T = x1 + (xC-x1)*lambda;

    y0T = y0 + dydx0*(xC-x0T);
    y1T = y1 + dydx1*(xC-x1T);

    dydxT0 = 0;
    dydxT1 = 0;

    yNegInf = 0;
    yInf    = sign(dydx1)*inf;    
    
    xScale = smoothness(j,1);

    if(j==1)
        yNegInf = yLimits(1,1);
        if(flag_appendCoeffcients==1)
            dydxT0 = dydx0;
            dydxT1 = dydx1;
        else
            dydxT0 = dydx0;
            dydxT1Err = dydx1-calcTanhSeriesDerivative(x1,...
                            tanhCoeffs,1);
            dydxT1 = curveInputs(j).dydx1+dydxT1Err;
            yInf   = sign(dydxT1)*inf;
        end
    else
        dydxT0 = 0;
        if(flag_appendCoeffcients==1)
            dydxT1 = dydx1-calcTanhSeriesDerivative(x1,...
                            tanhCoeffs(1:(j-1),:),1);
        else
            dydxT1Err = dydx1-calcTanhSeriesDerivative(x1,...
                            tanhCoeffs,1);
            dydxT1 = curveInputs(j).dydx1+dydxT1Err;
        end
        
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
        if(flag_appendCoeffcients==1)
            yPoint = y1 - calcTanhSeriesDerivative(xPoint,...
                            tanhCoeffs(1:(j-1),:),0);
        else
            yPointErr = y1 - calcTanhSeriesDerivative(xPoint,...
                            tanhCoeffs,0);
            yPoint =  curveInputs(j).yPoint + yPointErr;
        end
    end

    %disp(j);
    %fprintf('%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\n',...
    %    x0,x1,dydxT0,dydxT1,yNegInf,yInf,xScale,xPoint,yPoint);

    curveInputs(j).x0 = x0T;
    curveInputs(j).x1 = x1T;
    curveInputs(j).dydx0 = dydxT0;
    curveInputs(j).dydx1 = dydxT1;
    curveInputs(j).yNegInf= yNegInf;
    curveInputs(j).yInf = yInf;
    curveInputs(j).xScale = xScale;
    curveInputs(j).xPoint = xPoint;
    curveInputs(j).yPoint = yPoint;
    curveInputs(j).xAtIntYZero = xAtIntYZero;
    
    if(~((yPoint >= yNegInf && yPoint <= yInf && (yInf-yNegInf) > 0) ...
       || (yPoint <= yNegInf && yPoint >= yInf  && (yNegInf-yInf) > 0)) )
        here=1;
    end

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
