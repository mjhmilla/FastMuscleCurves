function figH = addDeGrooteFregly2016TendonCurveComparison(figH,...
                    curveParams,...
                    tendonForceLengthCurve, ...
                    ftDomainTest,...
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
expColor                = plotSettings.expColor;
expErrorColor           = plotSettings.expErrorColor;
tanColor                = plotSettings.tanColor;

eIso = curveParams.eIso;

c1 = 0.2;
c2 = 1.0;
c3 = 0.2;
m_kT = log((1.0 + c3) / c1) /(1.0 + eIso - c2);


npts = length(ftDomainTest);

ftN = zeros(npts,3);
ftBezierSample = zeros(npts,3);
ftErrorSample = zeros(npts,3);
ltN = ftDomainTest;

for i=1:1:3
    for j=1:1:length(ftDomainTest)
        
        ftN(j,i) = calcDeGrooteFregly2016TendonMultiplier(ltN(j,1),m_kT,c1,c2,c3,i-2);

        ftBezierSample(j,i) = calcBezierYFcnXDerivative(ltN(j,1),...
                                              tendonForceLengthCurve,i-2);
        ftErrorSample(j,i) = ftN(j,i)-ftBezierSample(j,i);        
    end

    figure(figH);
    subplot('Position',reshape(subPlotPanel(3,i,:),1,4));

    plot( ltN,ftBezierSample(:,i),...
          'Color',bezierColor,'LineWidth',1,...
          'DisplayName','Bezier');
    hold on;

    plot( ltN,ftN(:,i),...
          'Color',expColor,'LineWidth',1,...
          'DisplayName','DeGrooteFregly2016');
    hold on;

    yLimVec = ylim();
    dy = diff(yLimVec);

    text(ltN(end),...
         ftErrorSample(end,i)+dy*0.1,...
         sprintf('%1.2e: RMSE', sqrt(sum(ftErrorSample(:,i).^2))),...
         'FontSize',6,...
         'Color',expColor,...
         'HorizontalAlignment','right',...
         'VerticalAlignment','bottom');

    plot( ltN,ftErrorSample(:,i),...
          'Color',expErrorColor,'LineWidth',1,...
          'DisplayName','Tanh');
    hold on;


end


