function figH = addTanhTendonCurveComparison(figH,...
                        tendonCurveParams,...
                        tendonForceLengthCurve,...
                        ftDomainTest, ...
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
tanColor                = plotSettings.tanColor;

%%
% Create a C2 Quintic tendon force length curve
%%


x0      = tendonCurveParams.ltZeroN;
x1      = tendonCurveParams.ltToeN;
y0      = 0;
y1      = tendonCurveParams.fToeN;
dydx0   = 0;
dydx1   = tendonCurveParams.kToeN;




tanhSeriesParams(1).x0             = tendonCurveParams.ltZeroN;
tanhSeriesParams(1).x1             = tendonCurveParams.ltToeN;
tanhSeriesParams(1).dydx0          = 0;
tanhSeriesParams(1).dydx1          = tendonCurveParams.kToeN;
tanhSeriesParams(1).yNegInf        = 0;
tanhSeriesParams(1).yInf           = inf;
tanhSeriesParams(1).xScale         = 0.9;
tanhSeriesParams(1).xPoint         = tendonCurveParams.ltIsoN;
tanhSeriesParams(1).yPoint         = 1;
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

tanhSeriesCoefficients = [A,B,C,D,E,F];

optParams.names={'x0','x1','xScale'};
args        = [tanhSeriesParams(1).x0;...
               tanhSeriesParams(1).x1;...
               tanhSeriesParams(1).xScale];

argsScaling = args;
argsScaled  = ones(size(args));

errFcn = @(argInput)calcTanhCurveError(argInput,...
                 optParams,tanhSeriesParams,...
                 tendonForceLengthCurve,...
                 ftDomainTest,...
                 argsScaling);

errVec0 = errFcn(argsScaled);

[argScaledUpd,resnorm,residual,exitflag,output]=...
    lsqnonlin(errFcn,argsScaled);
argUpd = argScaledUpd.*argsScaling;

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

% yNegInf = y0;
% yInf    = [];  
% xShift      = -0.0045;
% xScale      = 0.125;
% xAtIntYZero = 1;
% 
% [A,B,C,D,E,F] = calcTanSegmentCoefficients(x0,x1,dydx0,dydx1,...
%                                            yNegInf,yInf, ...
%                                            xAtIntYZero,xShift,xScale);
% tendonForceLengthTanCoeffs = [A,B,C,D,E,F];

disp('The passive curves used in Brown, Scott, Loeb 1996 for the force-length');
disp('curves of the CE and the tendon are similar to the ones used here.');
disp('They should be referenced.');



npts=length(ftDomainTest);
ltN     = ftDomainTest;

tendonBezierSample      = zeros(npts,3);
tendonTanhSample        = zeros(npts,3);
tendonTanhError         = zeros(npts,3);
%tendonTanSample         = zeros(npts,3);


%%
%Plot the tendon curves
%%

for i=1:1:3
    for j=1:1:npts
        
        tendonBezierSample(j,i) = calcBezierYFcnXDerivative(ltN(j,1),...
                                              tendonForceLengthCurve,i-2);

        tendonTanhSample(j,i) = calcTanhSeriesDerivative(ltN(j,1),...
                                       tanhSeriesCoefficients,i-2);

%         tendonTanSample(j,i) = calcTanSeriesDerivative(ltN(j,1),...
%                                        tendonForceLengthTanCoeffs,i-2);

        tendonTanhError(j,i) = tendonTanhSample(j,i) ...
                             - tendonBezierSample(j,i); 

    end

    figure(figH);
    subplot('Position',reshape(subPlotPanel(indexPlotRow,i,:),1,4));

    if(flag_plotBezierCurves==1)
%         plot( ltN,tendonBezierSample(:,i),...
%               'Color',bezierColor,'LineWidth',bezierWidth,...
%               'DisplayName','Bezier');
%         hold on;
        fill([ltN(1,1);ltN(end,1);fliplr(ltN')'],...
             [0;0;fliplr(tendonBezierSample(:,i)')'],...
             bezierColor,...
             'EdgeColor','none',...
             'DisplayName','Bezier');
        hold on;
    end
    if(flag_plotTanhCurves==1)
        plot( ltN,tendonTanhSample(:,i),...
              'Color',tanhColor,'LineWidth',1,...
              'DisplayName','Tanh');
        hold on;
    end
%     if(flag_plotTanCurves==1)
%         plot( ltN,tendonTanSample(:,i),...
%               'Color',tanColor,'LineWidth',1,...
%               'DisplayName','Tan');
%         hold on;
%     end
    if(flag_plotTanhCurves==1 && flag_plotBezierCurves==1)
        text(ltN(end),...
             tendonTanhError(end,i),...
             sprintf('%1.2e: RMSE', sqrt(sum(tendonTanhError(:,i).^2))),...
             'FontSize',6,...
             'Color',tanhColor,...
             'HorizontalAlignment','right',...
             'VerticalAlignment','bottom');

        plot( ltN,tendonTanhError(:,i),...
              'Color',tanhErrorColor,'LineWidth',1,...
              'DisplayName','Tanh Error');
        hold on;
    end
    box off;

    switch(i)
        case 1
            title('Tendon Curve Integral');
            xlabel('Norm. Length ($\ell/\ell^T_s$)');
            ylabel('Norm. Energy ($\tilde{f}^T \, \tilde{\ell}^T$)');
            legend('Location','NorthEast');            
            legend boxoff;
        case 2
            plot([x0;x1],[y0;y1],'o','MarkerSize',5,...
                'Color',bezierColor,'MarkerFaceColor',[1,1,1]);
            hold on;
            title('Tendon Curve Value');            
            xlabel('Norm. Length ($\ell/\ell^T_s$)')
            ylabel('Norm. Force ($f/f^M_o$)');            

            
        case 3
            plot([x0;x1],[dydx0;dydx1],'o','MarkerSize',5,...
                'Color',bezierColor,'MarkerFaceColor',[1,1,1]);
            hold on;

            title('Derivative');                        
            xlabel('Norm. Length ($\ell/\ell^T_s$)')
            ylabel('Norm. Stiffness ($\delta \tilde{f^T} / \delta \tilde{\ell^T}$)');
            hold on;
            
        otherwise
            assert(0,'Error: missing postprocessing code for the current derivative');
    end

end