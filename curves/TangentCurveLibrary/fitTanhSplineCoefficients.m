function [tanhCoeffs,finalValuePolishing,finalDerivativePolishing] = ...
    fitTanhSplineCoefficients(xk,yk,dydxk,smoothness,...
    yLimits,dydxLimits,yLimitsSignOfAcceptableError,...
    xAtIntYZero,flag_polishKnotPoints,...
    valuesToPolish,derivativesToPolish,...
    maxPolishIterations, pinvTolerance,verbose)

assert(size(xk,2)==1,'Error: xk, yk, dydxk must be n x 1');
assert(length(xk)==length(yk) && length(xk)==length(dydxk),...
    'Error: xk, yk, and dydxk must be the same length');
assert(size(yLimits,1)==1 && size(yLimits,2)==2,...
    'Error: yLimits must be 1x2' );
assert(size(dydxLimits,1)==1 && size(dydxLimits,2)==2,...
    'Error: dydxLimits must be 1x2' );

n = length(xk)-1;


xS = [];
if(length(smoothness)==1)
    xS = ones(n,1).*smoothness;
else
    xS = smoothness;
end
dy = diff(yk).*0.05;

tanhCoeffs = zeros(n, 6);

%
% Construct the curve that will be close to passing through
% the knot points: xk,yk,dydxk
%

curveInputs(n) = struct('x0',0,'x1',0,'dydx0',0,'dydx1',0,...
                       'yNegInf',0,'yInf',0,'xScale',0,...
                       'xPoint',0,'yPoint',0,'xAtIntYZero',0);

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
    
    lambda = 0;
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

    flag_debug = 1;
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
        if(j < n)
            clf(figDebug);
        end

    end
end

%
% Polish up the knot points by making small adjustments to
%
%  B parameters to pass through (xk,yk)
%  D parameters to pass through xk, dydxk
%

xkVTarget = valuesToPolish(:,1);
xkDTarget = derivativesToPolish(:,1);

if(flag_polishKnotPoints==1)
    
      assert(isinf(valuesToPolish(1,1)) && isinf(valuesToPolish(end,1)),...
          ['Error: valuesToPolish must have x values of -inf and inf',...
          ' at the first and last indices']);

      assert(isinf(derivativesToPolish(1,1)) && isinf(derivativesToPolish(end,1)),...
          ['Error: derivativesToPolish must have x values of -inf and inf',...
          ' at the first and last indices']);

%     assert(indicesOfValueKnotsToPolish(1,1)==1 || isinf(yLimits(1,1)),...
%         ['Error: first and last y and dydx values must',...
%         ' be the first and last points to polish']);
% 
%     assert(indicesOfValueKnotsToPolish(1,end)==length(xk) || isinf(yLimits(1,2)),...
%         ['Error: first and last y and dydx values must',...
%         ' be the first and last points to polish']);
%     
%     assert(indicesOfDerivativeKnotsToPolish(1,1)==1,...
%         ['Error: first and last y and dydx values must',...
%         ' be the first and last points to polish']);
% 
%     assert(indicesOfDerivativeKnotsToPolish(1,end)==length(dydxk),...
%         ['Error: first and last y and dydx values must',...
%         ' be the first and last points to polish']);

    tanhCoeffsUpd = tanhCoeffs;

    % Approach
    %   Directly adjust B to fit yk, and A to fit dydxk
    arg = [tanhCoeffsUpd(:,1);tanhCoeffsUpd(:,2);tanhCoeffsUpd(:,3)];

    err = inf;
    iter =1;
    
    tol = eps*10*length(arg);
    iterMax=maxPolishIterations;
    errStart=0;


    stepLength=0.75;
    stepLengthMin = 0.01;
    stepLengthMax = 1;

    errBest=inf;
    tanhCoeffsBest=tanhCoeffs;
    argBest=arg;
    dargBest = ones(size(arg));
    darg = ones(size(arg)).*inf;

    xkVTarget       = valuesToPolish(:,1);%xk(indicesOfValueKnotsToPolish,1);
    ykVTarget       = valuesToPolish(:,2);%yk(indicesOfValueKnotsToPolish,1);    
    xkDTarget       = derivativesToPolish(:,1);%xk(indicesOfDerivativeKnotsToPolish,1);
    dydxkDTarget    = derivativesToPolish(:,2);%dydxk(indicesOfDerivativeKnotsToPolish,1);

    %Update the endpoints so that the values will be very close to the
    %numerical limit. Here I'm doing this in a lazy way by artifically
    %extending the minimum and maximum x points. Since the basis functions
    %are exponentials they will quickly converge to the limit.


    %Set the limits for the boundary evaluation inside the 
    %numerical limit when the regular function has to be exchanged for
    %the approximation.
    argLimit = log(realmax)*0.25;
    xkLim = [0,0];
    for k=1:1:size(tanhCoeffs,1)
        B = tanhCoeffsUpd(k,2);
        C = tanhCoeffsUpd(k,3);
        xkMin =-argLimit*C+B;
        xkMax = argLimit*C+B;

        if(xkMin < xkLim(1,1) || k==1)
            xkLim(1,1)=xkMin;
        end
        if(xkMax > xkLim(1,2) || k==1)
            xkLim(1,2) = xkMax;
        end
    end

    yLimitsSignOfAcceptableError;
    indexYEndpoints =[];
    signOfYEndpoints = [];

    for k=1:1:2
        if(~isinf(yLimits(1,k)))
            m=1;
            if(k==2)
                m = length(xkVTarget);
            end
            xkVTarget(m,1)      = xkLim(1,k);            
            ykVTarget(m,1)      = yLimits(1,k);
            indexYEndpoints     = [indexYEndpoints,m];
            if(m==1)
                signOfYEndpoints = [signOfYEndpoints,yLimitsSignOfAcceptableError(1,1)];
            else
                signOfYEndpoints = [signOfYEndpoints,yLimitsSignOfAcceptableError(1,2)];
            end
        else
            if(k==1)
                idxValid = [2:length(xkVTarget)]';
                xkVTarget = xkVTarget(idxValid,1);
                ykVTarget = ykVTarget(idxValid,1);                
            else
                idxValid = [1:(length(xkVTarget)-1)]';
                xkVTarget = xkVTarget(idxValid,1);
                ykVTarget = ykVTarget(idxValid,1);
            end
        end
        m=1;
        if(k==2)
            m = length(xkDTarget);
        end
        xkDTarget(m,1)      = xkLim(1,k);
        dydxkDTarget(m,1)   = dydxLimits(1,k);

    end
        

    
    nV =length(xkVTarget);
    nD = length(xkDTarget);
    m =length(tanhCoeffsUpd(:,2));

    errVM=inf;
    maxDivergingIterations = 20;
    iterDiverging = 0;

    while(errVM > tol && iter < iterMax && iterDiverging < maxDivergingIterations)
    
        %Update the tanh spline coefficients
        for i=1:1:m
            indexParams = [1;2;3];
            valueParams = [arg(i,1);arg(i+m,1);arg(i+2*m,1)];
    
            [A,B,C,D,E,F] = calcTanhSegmentCoefficientsGivenParameters(...
                curveInputs(i).x0,...
                curveInputs(i).x1,...
                curveInputs(i).dydx0,...
                curveInputs(i).dydx1,...
                curveInputs(i).yNegInf,...
                curveInputs(i).yInf, ...
                curveInputs(i).xScale,...
                curveInputs(i).xPoint, ...
                curveInputs(i).yPoint,...
                curveInputs(i).xAtIntYZero,...
                indexParams, ...
                valueParams);
            tanhCoeffsUpd(i,:) = [A,B,C,D,E,F];
        end

    
        ykC = zeros(size(xkVTarget));
        dydxkC = zeros(size(xkDTarget));
        
        for i=1:1:length(xkVTarget)
            ykC(i,1)   = calcTanhSeriesDerivative(...
                            xkVTarget(i,1),tanhCoeffsUpd,0);
        end
        for i=1:1:length(xkDTarget)
            dydxkC(i,1)= calcTanhSeriesDerivative(...
                            xkDTarget(i,1),tanhCoeffsUpd,1); 
        end    
           
        jacM = zeros(length(xkVTarget)+length(dydxkDTarget),length(arg));
        for i=1:1:length(xkVTarget)
            for j=1:1:m
                if(i==1 && j==3)
                    here=1;
                end

                jacM(i,j) = calcTanhSeriesParameterDerivative(...
                                xkVTarget(i,1),tanhCoeffsUpd(j,:),0,1);

                jacM(i,m+j) = calcTanhSeriesParameterDerivative(...
                                xkVTarget(i,1),tanhCoeffsUpd(j,:),0,2);

                jacM(i,2*m+j) = calcTanhSeriesParameterDerivative(...
                                xkVTarget(i,1),tanhCoeffsUpd(j,:),0,3);
               
            end
        end
        for i=1:1:length(dydxkDTarget)
            for j=1:1:m

                jacM(i+nV,j) = calcTanhSeriesParameterDerivative(...
                                xkDTarget(i,1),tanhCoeffsUpd(j,:),1,1);

                jacM(i+nV,m+j) = calcTanhSeriesParameterDerivative(...
                                xkDTarget(i,1),tanhCoeffsUpd(j,:),1,2);

                jacM(i+nV,2*m+j) = calcTanhSeriesParameterDerivative(...
                                xkDTarget(i,1),tanhCoeffsUpd(j,:),1,3);

            end
        end


        errV = [ykC;dydxkC] - [ykVTarget;dydxkDTarget];
        errVM = norm(errV);


        if(iter==1)
            errStart=errVM;
        end

        signsAcceptable=1;
        for z=1:1:length(indexYEndpoints)
            yEnd = ykC(indexYEndpoints(1,z),1);
            if(sign(yEnd)*signOfYEndpoints(1,z) < 0)
                signsAcceptable=0;
            end
        end

        if(errVM < errBest)
            if(signsAcceptable==1)
                errBest=errVM;
                tanhCoeffsBest=tanhCoeffsUpd;
                argBest=arg;
                dargBest=darg;
            end
            iterDiverging = 0;   

            jacMInv = pinv(jacM'*jacM,pinvTolerance);
            darg = -jacMInv*(jacM'*errV);
            arg = arg + darg*stepLength;              
        else
            iterDiverging=iterDiverging+1;            
        end



        %    stepLength=stepLength*1.2;
        %    stepLength=min(stepLength,stepLengthMax);
        %else
        %    stepLength=stepLength*0.5;
        %    stepLength = max(stepLength,stepLengthMin);
        %    arg= argBest + dargBest*stepLength;

        %end

        % f(q)+J(q)dq - f =0
        % J(q)dq = (f-f(q)
        % dq = (J(q)'*J(q)) \ (J(q)'*(f-f(q)))
        
        if(verbose==1)
            fprintf('%d\t%1.3e\t%1.3e\n',iter, errVM, stepLength);
        end


        iter=iter+1;
    end

    if(errBest > errStart)
        tanhCoeffsBest=tanhCoeffs;        
    else
        tanhCoeffs=tanhCoeffsBest;
    end

    here=1;
end


finalValuePolishing = zeros(length(xkVTarget),2);
finalDerivativePolishing = zeros(length(xkDTarget),2);

finalValuePolishing(:,1)=xkVTarget;
finalDerivativePolishing(:,1)=xkDTarget;

for i=1:1:length(xkVTarget)
    finalValuePolishing(i,2)  ...
        = calcTanhSeriesDerivative(finalValuePolishing(i,1),tanhCoeffs,0);
end
for i=1:1:length(xkDTarget)
    finalDerivativePolishing(i,2)  ...
        = calcTanhSeriesDerivative(finalDerivativePolishing(i,1),tanhCoeffs,1);
end  





