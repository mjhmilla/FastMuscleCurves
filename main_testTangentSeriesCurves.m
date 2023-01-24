clc;
close all;
clear all;

flag_plotBezierCurves=1;
flag_plotTanhCurves  =1;
flag_plotTanCurves   =0;

set(groot, 'defaultAxesFontSize',8);
set(groot, 'defaultTextFontSize',8);
set(groot, 'defaultAxesLabelFontSizeMultiplier',1.2);
set(groot, 'defaultAxesTitleFontSizeMultiplier',1.2);
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultAxesTitleFontWeight','bold');  
set(groot, 'defaultFigurePaperUnits','centimeters');
set(groot,'defaultFigurePaperType','A4');



pubOutputFolder                         = 'output/plots/MuscleCurves/';
postprocessingDirectoryTree             = genpath('postprocessing');
addpath(postprocessingDirectoryTree   );
addpath('colornames/');
[names,rgb] = colornames('SVG','FireBrick');

parametersDirectoryTreeMTParams     = genpath('parameters');
parametersDirectoryTreeExperiments  = genpath('experiments');
parametersDirectoryTreeModels       = genpath('models');
parametersDirectoryTreeCurves       = genpath('curves');

addpath(parametersDirectoryTreeMTParams);
addpath(parametersDirectoryTreeExperiments);
addpath(parametersDirectoryTreeModels);
addpath(parametersDirectoryTreeCurves);

%%
% Plotting configuration
%%

figCurves = figure;
bezierColor = [0,0,0];
tanhColor   = [0,0,1];
tanColor   = [1,0,1];

numberOfVerticalPlotRows      = 1;
numberOfHorizontalPlotColumns = 3;

plotWidth           = 6;
plotHeight          = 6;
plotHorizMarginCm   = 2;
plotVertMarginCm    = 2;

pageWidth = numberOfHorizontalPlotColumns*(plotWidth+plotHorizMarginCm)...
            + 1*plotHorizMarginCm;

pageHeight= numberOfVerticalPlotRows*(plotHeight+plotVertMarginCm)...
            + 1*plotVertMarginCm;

flag_usingOctave    = 0;


plotConfigGeneric;

%%
% Create a C2 Quintic tendon force length curve
%%

flag_enableNumericallyNonZeroGradients  = 0;
smallNumericallyNonZeroNumber           = sqrt(sqrt(eps));

eZero           = 1.0;
eIso            = 0.049;
kIso            = 1.375/eIso;
fToe            = 2./3.;
curvinessTendon = 0.5;
computeIntegral = 1;
minimumSlope    = 0;
if(flag_enableNumericallyNonZeroGradients==1)
  minimumSlope = smallNumericallyNonZeroNumber/10.;
end


%%
% Create and evaluate the two tendon curves
%%
tendonForceLengthCurve = ...
  createTendonForceLengthCurve2021( eZero, eIso, kIso, ...
                                    fToe, curvinessTendon, ...
                                    computeIntegral, ...
                                    flag_enableNumericallyNonZeroGradients,...
                                    smallNumericallyNonZeroNumber,...
                                    'TendonForceLengthCurve',...
                                    flag_usingOctave);

x0      = tendonForceLengthCurve.xEnd(1,1);
x1      = tendonForceLengthCurve.xEnd(1,2);
y0      = tendonForceLengthCurve.yEnd(1,1);
y1      = tendonForceLengthCurve.yEnd(1,2);
dydx0   = tendonForceLengthCurve.dydxEnd(1,1);
dydx1   = tendonForceLengthCurve.dydxEnd(1,2);


[A,B,C,D,E,F] = calcTanhSegmentCoefficients(x0,x1,dydx0,dydx1,y0,[],-0.00525,0.9,.0);
tendonForceLengthTanhCoeffs = [A,B,C,D,E,F];

[A,B,C,D,E,F] = calcTanSegmentCoefficients(x0,x1,dydx0,dydx1,y0,[], 1.0,-0.0045,0.125);
tendonForceLengthTanCoeffs = [A,B,C,D,E,F];

npts = 100;

ltN                     = zeros(npts,3);
tendonBezierSample      = zeros(npts,3);
tendonTanhSample        = zeros(npts,3);
tendonTanSample         = zeros(npts,3);


lNStart = x0-(x1-x0)*0.5;
lNEnd   = calcBezierFcnXGivenY(1,tendonForceLengthCurve,x1);
lNDelta = (lNEnd-lNStart)/(npts-1);
ltN     = [lNStart:lNDelta:lNEnd]';

for i=1:1:3
    for j=1:1:npts
        tendonBezierSample(j,i) = calcBezierYFcnXDerivative(ltN(j,1),...
                                              tendonForceLengthCurve,i-2);

        tendonTanhSample(j,i) = calcTanhSeriesDerivative(ltN(j,1),...
                                       tendonForceLengthTanhCoeffs,i-2);

        tendonTanSample(j,i) = calcTanSeriesDerivative(ltN(j,1),...
                                       tendonForceLengthTanCoeffs,i-2);

    end

    figure(figCurves);
    subplot('Position',reshape(subPlotPanel(1,i,:),1,4));

    if(flag_plotBezierCurves==1)
        plot( ltN,tendonBezierSample(:,i),...
              'Color',bezierColor,'LineWidth',1,...
              'DisplayName','Bezier');
        hold on;
    end
    if(flag_plotTanhCurves==1)
        plot( ltN,tendonTanhSample(:,i),...
              'Color',tanhColor,'LineWidth',1,...
              'DisplayName','Tanh');
        hold on;
    end
    if(flag_plotTanCurves==1)
        plot( ltN,tendonTanSample(:,i),...
              'Color',tanColor,'LineWidth',1,...
              'DisplayName','Tan');
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





rmpath(parametersDirectoryTreeMTParams);
rmpath(parametersDirectoryTreeExperiments);
rmpath(parametersDirectoryTreeModels);
rmpath(parametersDirectoryTreeCurves);

