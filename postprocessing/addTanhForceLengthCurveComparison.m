function figH = addTanhForceLengthCurveComparison(figH,...
                    curveParams, ...
                    fiberForceLengthCurve, ...
                    fpeDomainTest,...
                    plotSettings,...
                    flag_usingOctave)

subPlotPanel            = plotSettings.subPlotPanel;
indexPlotRow            = plotSettings.indexPlotRow;
flag_plotBezierCurves   = plotSettings.flag_plotBezierCurves;
flag_plotTanhCurves     = plotSettings.flag_plotTanhCurves;
flag_plotTanCurves      = plotSettings.flag_plotTanCurves;
bezierColor             = plotSettings.bezierColor;
tanhColor               = plotSettings.tanhColor;
tanhErrorColor          = plotSettings.tanhErrorColor;
tanColor                = plotSettings.tanColor;



% Create and evaluate the passive-force-length curves





tanhSeriesParams(1).x0             = curveParams.lpeZeroN;
tanhSeriesParams(1).x1             = curveParams.lpeToeN;
tanhSeriesParams(1).dydx0          = 0;
tanhSeriesParams(1).dydx1          = curveParams.kToeN;
tanhSeriesParams(1).yNegInf        = 0;
tanhSeriesParams(1).yInf           = inf;
tanhSeriesParams(1).xScale         = 0.9;
tanhSeriesParams(1).xPoint         = curveParams.lpeToeN;
tanhSeriesParams(1).yPoint         = curveParams.fToeN;
tanhSeriesParams(1).xAtIntYZero    = 0;

[A,B,C,D,E,F] = calcTanhSegmentCoefficientsUpd( ...
                    tanhSeriesParams(1).x0,     ...
                    tanhSeriesParams(1).x1,     ...
                    tanhSeriesParams(1).dydx0,  ...
                    tanhSeriesParams(1).dydx1,  ...
                    tanhSeriesParams(1).yNegInf,...
                    tanhSeriesParams(1).yInf,   ...
                    tanhSeriesParams(1).xScale, ...
                    tanhSeriesParams(1).xPoint, ...
                    tanhSeriesParams(1).yPoint, ...
                    tanhSeriesParams(1).xAtIntYZero);


forceLengthTanhCoeffs = [A,B,C,D,E,F];  


optParams.names={'x0','x1','xScale'};
args        = [tanhSeriesParams(1).x0;...
               tanhSeriesParams(1).x1;...
               tanhSeriesParams(1).xScale];

argScaling  = 1000;
argsScaled  = args .* argScaling;

errFcn = @(argInput)calcTanhCurveError(argInput,...
                 optParams,tanhSeriesParams,...
                 fiberForceLengthCurve,...
                 fpeDomainTest,...
                 argScaling);

errVec0 = errFcn(argsScaled);

[argScaledUpd,resnorm,residual,exitflag,output]=...
    lsqnonlin(errFcn,argsScaled);
argUpd = argScaledUpd./argScaling;

errVec1 = errFcn(argScaledUpd);

fprintf('%1.2e\tStarting Error\n%1.2e\tEnding Error\n',...
         sqrt(sum(errVec0.^2)),sqrt(sum(errVec1.^2)));  
fprintf('%i\tExit flag\n',exitflag);

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

    tanhSeriesCoefficients(i,:) = [A,B,C,D,E,F];
end

assert(length(optParams)==1);

x0 = localParams(1).x0;
x1 = localParams(1).x1;

y0 = calcTanhSeriesDerivative(x0,forceLengthTanhCoeffs,0);
y1 = calcTanhSeriesDerivative(x1,forceLengthTanhCoeffs,0);

dydx0 = calcTanhSeriesDerivative(x0,forceLengthTanhCoeffs,1);
dydx1 = calcTanhSeriesDerivative(x1,forceLengthTanhCoeffs,1);


%npts = 100;
%lceN = [0:(normLengthToe+0.1)/(npts-1):(normLengthToe+0.1)]';

lceN = fpeDomainTest;
npts = length(fpeDomainTest);

fpeBezierSample = zeros(npts,3);
fpeTanhSample   = zeros(npts,3);
fpeErrorSample = zeros(npts,3);

for i=1:1:3
    for j=1:1:npts
        fpeBezierSample(j,i) = calcBezierYFcnXDerivative(lceN(j,1),...
                                              fiberForceLengthCurve,i-2);

        fpeTanhSample(j,i) = calcTanhSeriesDerivative(lceN(j,1),...
                                       forceLengthTanhCoeffs,i-2);
        fpeErrorSample(j,i) = fpeTanhSample(j,i)-fpeBezierSample(j,i);
    end

    figure(figH);
    subplot('Position',reshape(subPlotPanel(2,i,:),1,4));

    if(flag_plotBezierCurves==1)
        plot( lceN,fpeBezierSample(:,i),...
              'Color',bezierColor,'LineWidth',1,...
              'DisplayName','Bezier');
        hold on;
    end
    if(flag_plotTanhCurves==1)
        plot( lceN,fpeTanhSample(:,i),...
              'Color',tanhColor,'LineWidth',1,...
              'DisplayName','Tanh');
        hold on;
    end
    if(flag_plotTanhCurves==1 && flag_plotBezierCurves==1)
        xp = lceN(end);        
        yp = fpeErrorSample(end,i);

        text( xp, yp,...
              sprintf('%1.2e: RMSE\n', sqrt(sum(fpeErrorSample(:,i).^2))),...
              'FontSize',6,...
              'HorizontalAlignment','right',...
              'VerticalAlignment','bottom',...
              'Color',tanhColor);
        hold on;
        plot( lceN,fpeErrorSample(:,i),...
              '-','Color',tanhErrorColor,'LineWidth',1,...
              'DisplayName','Tanh Error');
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
            title('Force-Length Curve Integral');
            xlabel('Norm. Length ($\ell/\ell^M_o$)');
            ylabel('Norm. Energy ($\tilde{f}^M \, \tilde{\ell}^M$)');
            legend('Location','NorthWest');            
            legend boxoff;
        case 2
            plot([x0;x1],[y0;y1],'o','MarkerSize',5,...
                'Color',bezierColor,'MarkerFaceColor',[1,1,1]);
            hold on;
            title('Force-Length  Curve Value');            
            xlabel('Norm. Length ($\ell/\ell^M_o$)')
            ylabel('Norm. Force ($f/f^M_o$)');            

            
        case 3
            plot([x0;x1],[dydx0;dydx1],'o','MarkerSize',5,...
                'Color',bezierColor,'MarkerFaceColor',[1,1,1]);
            hold on;

            title('Derivative');                        
            xlabel('Norm. Length ($\ell/\ell^M_o$)')
            ylabel('Norm. Stiffness ($\delta \tilde{f^M} / \delta \tilde{\ell^M}$)');
            hold on;
            
        otherwise
            assert(0,'Error: missing postprocessing code for the current derivative');
    end

end
