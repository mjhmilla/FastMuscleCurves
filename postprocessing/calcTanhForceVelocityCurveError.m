function errVector = calcTanhForceVelocityCurveError(argScaled, argNames, params,...
                        bezierCurve, domain, argScaling)

arg = argScaled ./ argScaling;
abcdefParams=zeros(length(params),6);

localParams = params;

idx=1;

for i=1:1:length(localParams)

    varNames = argNames(i).names;
    for j=1:1:length(varNames)
        localParams(i).(varNames{j})=arg(idx,1);
        idx=idx+1;
    end

    x0          =localParams(i).x0;
    x1          =localParams(i).x1;
    dydx0       =localParams(i).dydx0;
    dydx1       =localParams(i).dydx1;
    yNegInf     =localParams(i).yNegInf;
    yInf        =localParams(i).yInf;
    xScale      =localParams(i).xScale;
    xPoint      =localParams(i).xPoint;
    yPoint      =localParams(i).yPoint;
    xAtIntYZero =localParams(i).xAtIntYZero;

%     if(i==2)
%         dydxB = calcBezierYFcnXDerivative(1,bezierCurve,1);
%         dydxT = calcTanhSeriesDerivative(1,abcdefParams(1,:),1);
%         dydx1 = dydxB-dydxT;
%     end
    [A,B,C,D,E,F] = calcTanhSegmentCoefficientsUpd(...
                    x0,x1,dydx0,dydx1,...
                    yNegInf,yInf,...
                    xScale,xPoint, yPoint, xAtIntYZero);

    abcdefParams(i,:) = [A,B,C,D,E,F];
end

npts=100;

errVec= zeros(npts,1);
errDerVec= zeros(npts,1);

x0 = domain(1,1);
x1 = domain(1,2);
dx = (x1-x0)/(npts-1);

for i=1:1:npts
    xVal = x0 + dx*(i-1);
    val = calcBezierYFcnXDerivative(xVal,...
                                   bezierCurve,0);
    tanhVal = calcTanhSeriesDerivative(xVal,...
              abcdefParams,0);
    errVec(i,1)=tanhVal-val;

    dval = calcBezierYFcnXDerivative(xVal,...
                                   bezierCurve,1);
    dtanhVal = calcTanhSeriesDerivative(xVal,...
              abcdefParams,1);
    errDerVec(i,1)=dtanhVal-dval;
end

errVector = errVec;

flag_debug=0;
if(flag_debug==1)
    xVal        = [x0:dx:x1]';    
    yVal        = zeros(size(xVal));
    yValTanh    = zeros(size(xVal));

    for i=1:1:npts
        yVal(i,1) = calcBezierYFcnXDerivative(xVal(i,1),...
                                       bezierCurve,0);
        yValTanh(i,1) = calcTanhSeriesDerivative(xVal(i,1),...
                                abcdefParams,0);        
    end

    fig=figure;
    plot(xVal,yVal,'-k','DisplayName','Curve');
    hold on;
    plot(xVal,yValTanh,'-b','DisplayName','Tanh Curve');
    hold on;
    plot(xVal,errVec,'-r','DisplayName','Error');
    hold on;
    plot(xVal,yValTanh-yVal,'--k','DisplayName','Error (check)');
    hold on;
    xlabel('X');
    ylabel('Y');
    here=1;

end