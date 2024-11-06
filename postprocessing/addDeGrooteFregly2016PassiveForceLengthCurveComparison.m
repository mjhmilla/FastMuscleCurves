function figH = addDeGrooteFregly2016PassiveForceLengthCurveComparison(...
                    figH,...
                    curveParams,...
                    fiberForceLengthCurve, ...
                    fpeDomainTest,...
                    plotSettings)


subPlotPanel            = plotSettings.subPlotPanel;
indexPlotRow            = plotSettings.indexPlotRow;
flag_plotBezierCurves   = plotSettings.flag_plotBezierCurves;
flag_plotTanhCurves     = plotSettings.flag_plotTanhCurves;
flag_plotTanCurves      = plotSettings.flag_plotTanCurves;
bezierColor             = plotSettings.bezierColor;
tanhColor               = plotSettings.tanhColor;
tanhErrorColor          = plotSettings.tanhErrorColor;
expColor                = plotSettings.expColor;
expErrorColor           = plotSettings.expErrorColor;
tanColor                = plotSettings.tanColor;

e0 = curveParams.lceToeN-1;

kPE     = 4.0;
lceNMin = 0.2;

npts = length(fpeDomainTest);

fpeN = zeros(npts,3);
fpeBezierSample = zeros(npts,3);
fpeErrorSample = zeros(npts,3);
lceN = fpeDomainTest;

for i=1:1:3
    for j=1:1:length(fpeDomainTest)
        
        fpeN(j,i) = calcDeGrooteFregly2016PassiveMultiplier(lceN(j,1),...
                        e0,kPE,lceNMin,i-2);

        fpeBezierSample(j,i) = calcBezierYFcnXDerivative(lceN(j,1),...
                                              fiberForceLengthCurve,i-2);
        fpeErrorSample(j,i) = fpeN(j,i)-fpeBezierSample(j,i);        
    end

    figure(figH);
    subplot('Position',reshape(subPlotPanel(indexPlotRow,i,:),1,4));

    plot( lceN,fpeBezierSample(:,i),...
          'Color',bezierColor,'LineWidth',1,...
          'DisplayName','Bezier');
    hold on;

    plot( lceN,fpeN(:,i),...
          'Color',expColor,'LineWidth',1,...
          'DisplayName','DeGrooteFregly2016');
    hold on;

    yLimVec = ylim();
    dy = diff(yLimVec);

    text(lceN(end),...
         fpeErrorSample(end,i)+dy*0.1,...
         sprintf('%1.2e: RMSE', sqrt(sum(fpeErrorSample(:,i).^2))),...
         'FontSize',6,...
         'Color',expColor,...
         'HorizontalAlignment','right',...
         'VerticalAlignment','bottom');

    plot( lceN,fpeErrorSample(:,i),...
          'Color',expErrorColor,'LineWidth',1,...
          'DisplayName','DeGrooteFregly2016 error');
    hold on;


end


