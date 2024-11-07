function figH = addDeGrooteFregly2016ForceVelocityCurveComparison(...
                    figH,...
                    curveParams,...
                    forceVelocityCurve, ...
                    fvDomainTest,...
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

d1 = -0.3211346127989808;
d2 = -8.149;
d3 = -0.374;
d4 = 0.8825327733249912;

npts = length(fvDomainTest);

fvN = zeros(npts,1);
fvBezierSample = zeros(npts,1);
fvErrorSample = zeros(npts,1);
vceN = fvDomainTest;

idxCE = find(vceN >= -1 & vceN <= 1);
idxC = find(vceN >= -1 & vceN <= 0);
idxE = find(vceN >= 0 & vceN <= 1);

for j=1:1:length(fvDomainTest)
    
    fvN(j,1) = calcDeGrooteFregly2016ForceVelocityMultiplier(vceN(j,1),...
                    d1,d2,d3,d4,0);

    fvBezierSample(j,1) = calcBezierYFcnXDerivative(vceN(j,1),...
                                          forceVelocityCurve,0);
    fvErrorSample(j,1) = fvN(j,1)-fvBezierSample(j,1);        
end

figure(figH);
subplot('Position',reshape(subPlotPanel(indexPlotRow,2,:),1,4));

if(flag_plotBezierCurves==1)
%     plot( vceN,fvBezierSample(:,1),...
%           'Color',bezierColor,'LineWidth',bezierWidth,...
%           'DisplayName','Bezier');
%     hold on;
        fill([vceN(1,1);vceN(end,1);fliplr(vceN')'],...
             [0;0;fliplr(fvBezierSample(:,1)')'],...
             bezierColor,...
             'EdgeColor','none', 'DisplayName','Bezier');
        hold on;
end

plot( vceN,fvN(:,1),...
      'Color',expColor,'LineWidth',1,...
      'DisplayName','DeGrooteFregly2016');
hold on;

yLimVec = ylim();
dy = diff(yLimVec);


text(0,...
     fvErrorSample(end,1)+dy*0.1,...
     sprintf('%1.2e: RMSE\n%1.2e: Conc. RMSE\n%1.2e: Ecc. RMSE\n',...
          sqrt(sum(fvErrorSample(idxCE,1).^2)),...
          sqrt(sum(fvErrorSample(idxC,1).^2)),...
          sqrt(sum(fvErrorSample(idxE,1).^2)) ),...
          'FontSize',6,...
          'Color',expColor);
hold on;

plot( vceN,fvErrorSample(:,1),...
      'Color',expErrorColor,'LineWidth',1,...
      'DisplayName','DeGrooteFregly2016 error');
hold on;


box off;

