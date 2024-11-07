function figH = addDeGrooteFregly2016ActiveForceLengthCurveComparison(...
                    figH,...
                    curveParams,...
                    activeForceLengthCurve, ...
                    falDomainTest,...
                    plotSettings)


subPlotPanel            = plotSettings.subPlotPanel;
indexPlotRow            = plotSettings.indexPlotRow;
flag_plotBezierCurves   = plotSettings.flag_plotBezierCurves;
flag_plotTanhCurves     = plotSettings.flag_plotTanhCurves;
flag_plotTanCurves      = plotSettings.flag_plotTanCurves;
bezierColor             = plotSettings.bezierColor;
bezierWidth             = plotSettings.bezierWidth;
tanhColor               = plotSettings.tanhColor;
tanhErrorColor          = plotSettings.tanhErrorColor;
expColor                = plotSettings.expColor;
expErrorColor           = plotSettings.expErrorColor;
tanColor                = plotSettings.tanColor;

scale = curveParams.scale;

npts = length(falDomainTest);

falN            = zeros(npts,3);
falBezierSample = zeros(npts,3);
falErrorSample  = zeros(npts,3);
lceN            = falDomainTest;

for i=2:1:3
    for j=1:1:length(falDomainTest)
        
        falN(j,i) = calcDeGrooteFregly2016ActiveForceLengthMultiplier(...
                        lceN(j,1),scale,i-2);

        falBezierSample(j,i) = calcBezierYFcnXDerivative(lceN(j,1),...
                                activeForceLengthCurve,i-2);
        falErrorSample(j,i) = falN(j,i)-falBezierSample(j,i);        
    end

    figure(figH);
    subplot('Position',reshape(subPlotPanel(indexPlotRow,i,:),1,4));

    if(flag_plotBezierCurves==1)
%         plot( lceN,falBezierSample(:,i),...
%               'Color',bezierColor,'LineWidth',bezierWidth,...
%               'DisplayName','Bezier');
%         hold on;
        fill([lceN(1,1);lceN(end,1);fliplr(lceN')'],...
             [0;0;fliplr(falBezierSample(:,i)')'],...
             bezierColor,...
             'EdgeColor','none');
        hold on;
    end

    plot( lceN,falN(:,i),...
          'Color',expColor,'LineWidth',1,...
          'DisplayName','DeGrooteFregly2016');
    hold on;

    yLimVec = ylim();
    dy = diff(yLimVec);

    text(max(lceN),...
         0.95*max(falN(:,i)),...
         sprintf('%1.2e: RMSE', sqrt(sum(falErrorSample(:,i).^2))),...
         'FontSize',6,...
         'Color',expColor,...
         'HorizontalAlignment','right',...
         'VerticalAlignment','top');

    plot( lceN,falErrorSample(:,i),...
          'Color',expErrorColor,'LineWidth',1,...
          'DisplayName','DeGrooteFregly2016 error');
    hold on;

    box off;

end


