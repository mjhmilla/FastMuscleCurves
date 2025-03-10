function tanhSeriesCurves = createTanhSeriesCurvesFromBezierCurves(...
                              bezierCurves, minYLimit, verbose)

bezierCurveNames = fields(bezierCurves);
pinvTolerance=1e-8;
maxPolishInterations = 500;

for indexCurve=1:1:length(bezierCurveNames)

    curveName = bezierCurveNames{indexCurve};

    %xk = [bezierCurves.(curveName).xpts(1,1:(end))';...
    %      bezierCurves.(curveName).xpts(end,end)];
    xk = bezierCurves.(curveName).xk;

    numberOfSubdivisions = bezierCurves.(curveName).numberOfSubdivisions;
    if(numberOfSubdivisions>0)
        xkUpd = zeros((numberOfSubdivisions+1)*(length(xk)-1)+1,1);
        smoothnessUpd = zeros(length(xkUpd)-1,1);
        idx=1;
        for i=1:1:(length(xk)-1)
            for j=1:1:(1+numberOfSubdivisions)
                x0 = xk(i,1);
                x1 = xk(i+1,1);
                delta = (j-1)/(numberOfSubdivisions+1);
                xkUpd(idx,1)=x0 + (x1-x0)*delta;

                if((i+1)<=length(bezierCurves.(curveName).smoothness))
                    s0 = bezierCurves.(curveName).smoothness(i,1);
                    s1 = bezierCurves.(curveName).smoothness(i+1,1);  
                    smoothnessUpd(idx,1) = s0 + (s1-s0)*delta;
                else
                    smoothnessUpd(idx,1) = smoothnessUpd(idx-1,1);
                end

                idx=idx+1;
            end
        end
        xkUpd(idx,1)=xk(end,1);
        xk=xkUpd;
        bezierCurves.(curveName).smoothness=smoothnessUpd;
    end
    yk      = zeros(size(xk));
    dydxk   = zeros(size(xk));
    for i=1:1:length(xk)
        yk(i,1)=calcBezierYFcnXDerivative(xk(i,1),bezierCurves.(curveName),0);
        dydxk(i,1)=calcBezierYFcnXDerivative(xk(i,1),bezierCurves.(curveName),1);
    end
    yk(1,1)=bezierCurves.(curveName).yEnd(1,1);
    yk(end,1)=bezierCurves.(curveName).yEnd(1,2);    
    dydxk(1,1)=bezierCurves.(curveName).dydxEnd(1,1);
    dydxk(end,1)=bezierCurves.(curveName).dydxEnd(1,2);

    yLimits    = bezierCurves.(curveName).yEnd;
    for z=1:1:length(yLimits)
        yLimits(1,z) = max(minYLimit,yLimits(1,z));
    end
    

    yLimitsSignOfAcceptableError =...
        bezierCurves.(curveName).yLimitsSignOfAcceptableError; 
    %In this case all curves must be above zero
    dydxLimits = bezierCurves.(curveName).dydxEnd;
    for z=1:1:length(dydxLimits)
        if(abs(dydxLimits(1,z))>eps)
            yLimits(1,z) = sign(dydxLimits(1,z))*inf;
        end
    end

    nSample=100;
    x0 = bezierCurves.(curveName).xEnd(1,1);
    x1 = bezierCurves.(curveName).xEnd(1,2);
    xSample = [x0:((x1-x0)/(nSample-1)):x1]';

    xSample = unique([xSample;xk]);
    xSample = sort(xSample);
    %     ySample=zeros(size(xSample));
%     for i=1:1:length(xSample)
%         ySample(i,1) = calcBezierYFcnXDerivative(xSample(i,1),...
%                          bezierCurves.activeForceLengthCurve,0);
%     end

    xAtIntYZero = bezierCurves.(curveName).xAtIntYZero;    
    flag_polishKnotPoints=1;
    %indexOfValueKnotsToPolish = ...
    %    bezierCurves.(curveName).indexOfValueKnotsToPolish;
    %indexOfDerivativeKnotsToPolish = ...
    %    bezierCurves.(curveName).indexOfDerivativeKnotsToPolish;

    [tanhSplineCoeffs, finalPolishedValues, finalPolishedDerivatives ]= ...
        fitTanhSplineCoefficients(...
            xk,yk,dydxk,bezierCurves.(curveName).smoothness,...
            yLimits,dydxLimits,yLimitsSignOfAcceptableError,...
            xAtIntYZero,flag_polishKnotPoints,...
            bezierCurves.(curveName).valuesToPolish,...
            bezierCurves.(curveName).firstDerivativeValuesToPolish,...
            maxPolishInterations,...
            pinvTolerance,...
            verbose);
    
    ySample = zeros(size(xSample));
    dySample = zeros(size(xSample));
    yIntSample = zeros(size(xSample));
    
    for i=1:1:length(xSample)
        ySample(i,1)= ...
            calcTanhSeriesDerivative(xSample(i,1),tanhSplineCoeffs,0);
        dySample(i,1)= ...
            calcTanhSeriesDerivative(xSample(i,1),tanhSplineCoeffs,1);
        yIntSample(i,1)= ...
            calcTanhSeriesDerivative(xSample(i,1),tanhSplineCoeffs,-1);
    end

    tanhSeriesCurves.(curveName).curveParameters=tanhSplineCoeffs;
    tanhSeriesCurves.(curveName).knots.x = xk;
    tanhSeriesCurves.(curveName).knots.y = yk;    
    tanhSeriesCurves.(curveName).sample.x = xSample;
    tanhSeriesCurves.(curveName).sample.y = ySample;
    tanhSeriesCurves.(curveName).sample.dy = dySample;
    tanhSeriesCurves.(curveName).sample.yInt = yIntSample;
    tanhSeriesCurves.(curveName).finalPolishedValues=finalPolishedValues;
    tanhSeriesCurves.(curveName).finalPolishedDerivatives=finalPolishedDerivatives;

    fprintf('\n\nPolished values and derivatives\n\n');
    for z=1:1:length(finalPolishedValues)
        fprintf('%i\t%1.3e\t%1.3e\tvalue\n',...
            z,...
            bezierCurves.(curveName).valuesToPolish(z,2),...
            finalPolishedValues(z,2));
        here=1;
    end
    for z=1:1:length(finalPolishedDerivatives)
        fprintf('%i\t%1.3e\t%1.3e\tdydx value\n',...
            z,...
            bezierCurves.(curveName).firstDerivativeValuesToPolish(z,2),...
            finalPolishedDerivatives(z,2));
        here=1;
    end

    here=1;
end