function figH = addTanhForceVelocityCurveComparison01(figH,...
    plotSettings, flag_usingOctave)

flag_fitfvCurve=1;
flag_plotComponents=0;
subPlotPanel            = plotSettings.subPlotPanel;
indexPlotRow            = plotSettings.indexPlotRow;
flag_plotBezierCurves   = plotSettings.flag_plotBezierCurves;
flag_plotTanhCurves     = plotSettings.flag_plotTanhCurves;
flag_plotTanCurves      = plotSettings.flag_plotTanCurves;
bezierColor             = plotSettings.bezierColor;
tanhColor               = plotSettings.tanhColor;
tanColor                = plotSettings.tanColor;


flag_enableNumericallyNonZeroGradients  = 0;
smallNumericallyNonZeroNumber           = sqrt(sqrt(eps));
smallNumericallyNonZeroSlope            = sqrt(eps);

curvinessEccentricForceVelocity = 1.0;
flag_sharpEccentricTransition = 0;

forceVelocityMultiplierAtHalfMaximumFiberVelocity       = 0.2;
forceVelocityMultiplierAtLowEccentricFiberVelocity      = 1.4;
forceVelocityMultiplierAtMaximumEccentricFiberVelocity  = 1.5;

fiberForceVelocityCurve ...
  = createFiberForceVelocityCurve(...
      forceVelocityMultiplierAtHalfMaximumFiberVelocity,...
      forceVelocityMultiplierAtLowEccentricFiberVelocity,...
      forceVelocityMultiplierAtMaximumEccentricFiberVelocity,...
      curvinessEccentricForceVelocity,...
      flag_sharpEccentricTransition,...
      flag_enableNumericallyNonZeroGradients,...
      smallNumericallyNonZeroNumber,...
      smallNumericallyNonZeroSlope,...
      'ForceVelocityCurve',...
      flag_usingOctave);

tanhSeriesParams(2) = ...
    struct('x0',0,'dydx0',0,...
           'x1',0,'dydx1',0,...
           'yNegInf',0,'yInf',0,...
           'xScale',0,'xPoint',0,...
           'yPoint',0,'xAtIntYZero',0);

optParams(2) = struct('names',{''});

dyPoint = 0.01;
yInf    = inf;
yNegInf = 0;
wC = 0.4;

idxT = 1;

tanhSeriesParams(idxT).x0          = -1.5;
tanhSeriesParams(idxT).x1          = 0;
tanhSeriesParams(idxT).dydx0       = 0;
tanhSeriesParams(idxT).dydx1       = ...
    0.25*wC*calcBezierYFcnXDerivative(  ...
        tanhSeriesParams(idxT).x1, ...
        fiberForceVelocityCurve,1);
tanhSeriesParams(idxT).yNegInf     = yNegInf;
tanhSeriesParams(idxT).yInf        = yInf;
tanhSeriesParams(idxT).xScale      = 1.0;
tanhSeriesParams(idxT).xPoint      = tanhSeriesParams(idxT).x1;
tanhSeriesParams(idxT).yPoint      = ...
      wC*(calcBezierYFcnXDerivative(  ...
        tanhSeriesParams(idxT).x1, ...
        fiberForceVelocityCurve,0) ...
        + dyPoint);
tanhSeriesParams(idxT).xAtIntYZero = tanhSeriesParams(idxT).x0;

optParams(idxT).names = {'x0','xScale','dydx1'};

args = [tanhSeriesParams(idxT).x0;...
        tanhSeriesParams(idxT).xScale; ...
        tanhSeriesParams(idxT).dydx1];

[A,B,C,D,E,F] = calcTanhSegmentCoefficientsUpd( ...
                    tanhSeriesParams(idxT).x0,     ...
                    tanhSeriesParams(idxT).x1,     ...
                    tanhSeriesParams(idxT).dydx0,  ...
                    tanhSeriesParams(idxT).dydx1,  ...
                    tanhSeriesParams(idxT).yNegInf,...
                    tanhSeriesParams(idxT).yInf,   ...
                    tanhSeriesParams(idxT).xScale, ...
                    tanhSeriesParams(idxT).xPoint, ...
                    tanhSeriesParams(idxT).yPoint, ...
                    tanhSeriesParams(idxT).xAtIntYZero);

forceVelocityTanhCoeffs = [A,B,C,D,E,F];  

tanhSeriesParams(2).x0          = -0.6;
tanhSeriesParams(2).x1          = 0;
tanhSeriesParams(2).dydx0       = 0;
tanhSeriesParams(2).dydx1       = ...
    (1-wC)*calcBezierYFcnXDerivative(  ...
            tanhSeriesParams(2).x1, ...
            fiberForceVelocityCurve,1);
tanhSeriesParams(2).yNegInf     = yNegInf;
tanhSeriesParams(2).yInf        = yInf;
tanhSeriesParams(2).xScale      = 1.0;
tanhSeriesParams(2).xPoint      = tanhSeriesParams(idxT).x1;
tanhSeriesParams(2).yPoint      = ...
    (1-wC)*(calcBezierYFcnXDerivative(  ...
            tanhSeriesParams(2).x1, ...
            fiberForceVelocityCurve,0) ...
            + dyPoint);
tanhSeriesParams(2).xAtIntYZero = tanhSeriesParams(1).x0;

optParams(2).names = {'x0','xScale','dydx1'};

args = [args;...
        tanhSeriesParams(2).x0;...
        tanhSeriesParams(2).xScale; ...
        tanhSeriesParams(2).dydx1];

[A,B,C,D,E,F] = calcTanhSegmentCoefficientsUpd( ...
                    tanhSeriesParams(2).x0,     ...
                    tanhSeriesParams(2).x1,     ...
                    tanhSeriesParams(2).dydx0,  ...
                    tanhSeriesParams(2).dydx1,  ...
                    tanhSeriesParams(2).yNegInf,...
                    tanhSeriesParams(2).yInf,   ...
                    tanhSeriesParams(2).xScale, ...
                    tanhSeriesParams(2).xPoint, ...
                    tanhSeriesParams(2).yPoint, ...
                    tanhSeriesParams(2).xAtIntYZero);

forceVelocityTanhCoeffs = [forceVelocityTanhCoeffs;...
                           A,B,C,D,E,F];  


%%
%Get the slope at postive infinity
%%
flag_addEccentric=1;
if(flag_addEccentric==1)
    %Add a curve that will flatten the eccentric side to the desired final
    %slope
    
    tanhSeriesParams(3).x0          = 0;
    tanhSeriesParams(3).x1          = 0.5;
    tanhSeriesParams(3).dydx0       = 0;
    tanhSeriesParams(3).dydx1       = ...
        (forceVelocityMultiplierAtMaximumEccentricFiberVelocity...
        -forceVelocityMultiplierAtLowEccentricFiberVelocity)/(0.9) ...
        - calcTanhSeriesDerivative(1e6,forceVelocityTanhCoeffs,1);
    tanhSeriesParams(3).yNegInf     = 0;
    tanhSeriesParams(3).yInf        = sign(tanhSeriesParams(3).dydx1)*yInf;
    tanhSeriesParams(3).xScale      = 1;
    tanhSeriesParams(3).xPoint      = 0;
    tanhSeriesParams(3).yPoint      = -dyPoint;
    tanhSeriesParams(3).xAtIntYZero = 0;
    
    optParams(3).names = {'x1','xScale'};
    args = [args; 
           tanhSeriesParams(3).x1;...
           tanhSeriesParams(3).xScale];
    
    [A,B,C,D,E,F] = calcTanhSegmentCoefficientsUpd( ...
                        tanhSeriesParams(3).x0,     ...
                        tanhSeriesParams(3).x1,     ...
                        tanhSeriesParams(3).dydx0,  ...
                        tanhSeriesParams(3).dydx1,  ...
                        tanhSeriesParams(3).yNegInf,...
                        tanhSeriesParams(3).yInf,   ...
                        tanhSeriesParams(3).xScale, ...
                        tanhSeriesParams(3).xPoint, ...
                        tanhSeriesParams(3).yPoint, ...
                        tanhSeriesParams(3).xAtIntYZero);
    
    forceVelocityTanhCoeffs = [forceVelocityTanhCoeffs;...
                               A,B,C,D,E,F];  
end

optDomain = [-1,1];
argScaling = 1000;
argsScaled = args.*argScaling;

errVec0 = calcTanhForceVelocityCurveError(argsScaled, optParams,...
            tanhSeriesParams,fiberForceVelocityCurve, optDomain,...
            argScaling);

errFcn = @(argInput)calcTanhForceVelocityCurveError(argInput,optParams,...
           tanhSeriesParams,fiberForceVelocityCurve,optDomain,argScaling);

if(flag_fitfvCurve==1)

    [argScaledUpd,resnorm,residual,exitflag,output]=lsqnonlin(errFcn,argsScaled);
    
    errVec1 = calcTanhForceVelocityCurveError(argScaledUpd, optParams,...
                tanhSeriesParams,fiberForceVelocityCurve, optDomain,argsScaled);

    argUpd = argScaledUpd ./ argScaling;

    fprintf('%1.2e\tStarting Error\n%1.2e\tEnding Error',...
             sqrt(sum(errVec0.^2)),sqrt(sum(errVec1.^2)));
    
    localParams=tanhSeriesParams;
    idx=1;
    for i=1:1:length(optParams)
    
        varNames = optParams(i).names;
        for j=1:1:length(varNames)
            localParams(i).(varNames{j})=argUpd(idx,1);
            idx=idx+1;
        end
    
        x0_          =localParams(i).x0;
        x1_          =localParams(i).x1;
        dydx0_       =localParams(i).dydx0;
        dydx1_       =localParams(i).dydx1;
        yNegInf_     =localParams(i).yNegInf;
        yInf_        =localParams(i).yInf;
        xScale_      =localParams(i).xScale;
        xPoint_      =localParams(i).xPoint;
        yPoint_      =localParams(i).yPoint;
        xAtIntYZero_ =localParams(i).xAtIntYZero;
    
        [A,B,C,D,E,F] = calcTanhSegmentCoefficientsUpd(...
                        x0_,x1_,dydx0_,dydx1_,...
                        yNegInf_,yInf_,...
                        xScale_,xPoint_, yPoint_, xAtIntYZero_);
    
        forceVelocityTanhCoeffs(i,:) = [A,B,C,D,E,F];
    end
end

here=1;

npts = 301;
xP0 = -1.5;
xP1 = 1.5;
vceN = [(xP0):(xP1-xP0)/(npts-1):(xP1)]';

idxCE = find(vceN >= -1 & vceN <= 1);
idxC = find(vceN >= -1 & vceN <= 0);
idxE = find(vceN >= 0 & vceN <= 1);

fvBezierSample = zeros(npts,3);
fvErrorSample = zeros(npts,3);

fvTanhSample        = zeros(npts,3);
fvComponentSample   = zeros(npts,size(forceVelocityTanhCoeffs,1));

for i=2:1:3
    for j=1:1:npts
        fvBezierSample(j,i) = calcBezierYFcnXDerivative(vceN(j,1),...
                                              fiberForceVelocityCurve,i-2);

        fvTanhSample(j,i) = calcTanhSeriesDerivative(vceN(j,1),...
                                       forceVelocityTanhCoeffs,i-2);

        fvErrorSample(j,i)=fvTanhSample(j,i)-fvBezierSample(j,i);

%         if(i==2)
%             y0 =0.1+0.1*tanh((0-0.05)/0.01);
%             fvTanhSample(j,i) = fvTanhSample(j,i)+ 0.1+0.1*tanh((vceN(j,1)-0.05)/0.01)-y0;
%         end

        if(i-2 ==0)
            for k=1:1:size(fvComponentSample,2)
                fvComponentSample(j,k) = calcTanhSeriesDerivative(vceN(j,1),...
                                       forceVelocityTanhCoeffs(k,:),i-2);
            end
        end
    end

    figure(figH);
    subplot('Position',reshape(subPlotPanel(indexPlotRow ,i,:),1,4));

    if(flag_plotBezierCurves==1)
        plot( vceN,fvBezierSample(:,i),...
              'Color',bezierColor,'LineWidth',1,...
              'DisplayName','Bezier');
        hold on;
    end
    if(flag_plotTanhCurves==1)
        plot( vceN,fvTanhSample(:,i),...
              'Color',tanhColor,'LineWidth',1,...
              'DisplayName','Tanh');
        hold on;
        if(i==2 && flag_plotComponents==1)
            for k=1:1:size(fvComponentSample,2)
                plot( vceN,fvComponentSample(:,k),...
                      'Color',tanhColor.*0.5+[1,1,1].*0.5,'LineWidth',1,...
                      'DisplayName','Tanh-Component');
                hold on;
            end
        end
    end
    if(flag_plotBezierCurves==1 && flag_plotTanhCurves==1)
        plot(vceN,fvErrorSample(:,i),'-r','DisplayName','Error');
        hold on;
        text(-1,1,sprintf('%1.2e: RMSE\n%1.2e: Conc. RMSE\n%1.2e: Ecc. RMSE\n',...
                  sqrt(sum(fvErrorSample(idxCE,i).^2)),...
                  sqrt(sum(fvErrorSample(idxC,i).^2)),...
                  sqrt(sum(fvErrorSample(idxE,i).^2)) ),...
                  'FontSize',6);
        hold on;

    end
%     if(flag_plotTanCurves==1)
%         plot( lceN,fpeTanSample(:,i),...
%               'Color',tanColor,'LineWidth',1,...
%               'DisplayName','Tan');
%         hold on;
%     end
    box off;

    switch(i)
        case 1
            title('Force-Velocity Curve Integral');
            xlabel('Norm. Velocity ($v/v^M_{max}$)');
            ylabel('Norm. Int. Force-Velocity ($\int \tilde{f}^M \tilde{v}^M$)');
            legend('Location','NorthEast');            
            legend boxoff;
        case 2
            plot([tanhSeriesParams(1).x0;...
                  tanhSeriesParams(1).x1],...
                  [0;...
                   tanhSeriesParams(1).yPoint],'o','MarkerSize',5,...
                'Color',bezierColor,'MarkerFaceColor',[1,1,1]);
            hold on;
            title('Force-Velocity Curve Value');            
            xlabel('Norm. Velocity ($v/v^M_{max}$)')
            ylabel('Norm. Force ($f/f^M_o$)');            

            
        case 3
            plot([tanhSeriesParams(1).x0;...
                  tanhSeriesParams(1).x1],...
                  [tanhSeriesParams(1).dydx0;...
                   tanhSeriesParams(1).dydx1],'o','MarkerSize',5,...
                'Color',bezierColor,'MarkerFaceColor',[1,1,1]);
            hold on;

            title('Derivative');                        
            xlabel('Norm. Velocity ($v/v^M_{max}$)')
            ylabel('Norm. Slope ($(f/f^M_o)/(v/v^M_{max})$)');
            hold on;
            
        otherwise
            assert(0,'Error: missing postprocessing code for the current derivative');
    end

end