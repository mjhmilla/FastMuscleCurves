clc;
close all;
clear all;

flag_plotBezierCurves=1;
flag_plotTanhCurves  =1;
flag_plotTanCurves   =0;

flag_enableNumericallyNonZeroGradients=1;
smallNumericallyNonZeroValue=sqrt(eps);
smallNumericallyNonZeroSlope=sqrt(sqrt(eps));

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

numberOfVerticalPlotRows      = 3;
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


[subPlotPanel,pageWidth,pageHeight]  = ...
    plotConfigGeneric(numberOfHorizontalPlotColumns, ...
                      numberOfVerticalPlotRows,...
                      plotWidth,...
                      plotHeight,...
                      plotHorizMarginCm,...
                      plotVertMarginCm);
%^plotConfigGeneric;

%%
% Create a C2 Quintic tendon force length curve
%%

flag_enableNumericallyNonZeroGradients  = 0;
smallNumericallyNonZeroNumber           = sqrt(sqrt(eps));
smallNumericallyNonZeroSlope            = sqrt(eps);
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
                                    smallNumericallyNonZeroSlope,...
                                    'TendonForceLengthCurve',...
                                    flag_usingOctave);

x0      = tendonForceLengthCurve.xEnd(1,1);
x1      = tendonForceLengthCurve.xEnd(1,2);
y0      = tendonForceLengthCurve.yEnd(1,1);
y1      = tendonForceLengthCurve.yEnd(1,2);
dydx0   = tendonForceLengthCurve.dydxEnd(1,1);
dydx1   = tendonForceLengthCurve.dydxEnd(1,2);

yNegInf = y0;
yInf    = inf;  

xShift      = -0.00525;
xScale      = 0.9;
xAtIntYZero = 0;
xPoint = 1+eIso;
yPoint = 1;

[A,B,C,D,E,F] = calcTanhSegmentCoefficientsUpd(...
                    x0,x1,dydx0,dydx1,...
                    yNegInf,yInf,...
                    xScale,xPoint, yPoint, xAtIntYZero);

%[A,B,C,D,E,F] = calcTanhSegmentCoefficients(x0,x1,dydx0,dydx1,...
%                                            yNegInf,yInf,...
%                                            xShift,xScale,xAtIntYZero);
tendonForceLengthTanhCoeffs = [A,B,C,D,E,F];

yNegInf = y0;
yInf    = [];  
xShift      = -0.0045;
xScale      = 0.125;
xAtIntYZero = 1;



[A,B,C,D,E,F] = calcTanSegmentCoefficients(x0,x1,dydx0,dydx1,...
                                           yNegInf,yInf, ...
                                           xAtIntYZero,xShift,xScale);
tendonForceLengthTanCoeffs = [A,B,C,D,E,F];

disp('The passive curves used in Brown, Scott, Loeb 1996 for the force-length');
disp('curves of the CE and the tendon are similar to the ones used here.');
disp('They should be referenced.');
tendonTanhToe = calcTanhSeriesDerivative(x1,tendonForceLengthTanhCoeffs,0);
errToe = tendonTanhToe;


npts = 100;

ltN                     = zeros(npts,3);
tendonBezierSample      = zeros(npts,3);
tendonTanhSample        = zeros(npts,3);
tendonTanSample         = zeros(npts,3);


lNStart = x0-(x1-x0)*0.5;
lNEnd   = calcBezierFcnXGivenY(1,tendonForceLengthCurve,x1);
lNDelta = (lNEnd-lNStart)/(npts-1);
ltN     = [lNStart:lNDelta:lNEnd]';

%%
%Plot the tendon curves
%%

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

%%
% Create and evaluate the passive-force-length curves
%%
  normLengthZero = 1+0; 
  normLengthToe  = 1+0.6;
  fToe  = 1;
  yZero = 0;
  kZero = 0;
  if(flag_enableNumericallyNonZeroGradients)
    yZero   = smallNumericallyNonZeroValue;
    kZero   = smallNumericallyNonZeroSlope; 
  end         
  kLow  = 0.2;
  kToe  = 2/(normLengthToe-normLengthZero);
  curviness = 0.75;  
  flag_computeIntegral = 1;

  fiberForceLengthCurve = ...
    createFiberForceLengthCurve2021(normLengthZero,...
                                normLengthToe,...
                                fToe,...
                                yZero,...
                                kZero,...
                                kLow,...
                                kToe,...
                                curviness,...
                                flag_computeIntegral,...
                                'FiberForceLengthCurve',...
                                flag_usingOctave);   


x0      = fiberForceLengthCurve.xEnd(1,1);
x1      = fiberForceLengthCurve.xEnd(1,2);
y0      = fiberForceLengthCurve.yEnd(1,1);
y1      = fiberForceLengthCurve.yEnd(1,2);
dydx0   = fiberForceLengthCurve.dydxEnd(1,1);
dydx1   = fiberForceLengthCurve.dydxEnd(1,2);

yNegInf = y0;
yInf    = inf;  

xScale      = 0.9;
xAtIntYZero = 0;
xPoint = normLengthToe;
yPoint = fToe;

[A,B,C,D,E,F] = calcTanhSegmentCoefficientsUpd(...
                    x0,x1,dydx0,dydx1,...
                    yNegInf,yInf,...
                    xScale,xPoint, yPoint, xAtIntYZero);

%[A,B,C,D,E,F] = calcTanhSegmentCoefficients(x0,x1,dydx0,dydx1,...
%                                            yNegInf,yInf,...
%                                            xShift,xScale,xAtIntYZero);
forceLengthTanhCoeffs = [A,B,C,D,E,F];  

npts = 100;
lceN = [0:(normLengthToe+0.1)/(npts-1):(normLengthToe+0.1)]';

fpeBezierSample = zeros(npts,3);
fpeTanhSample   = zeros(npts,3);

for i=1:1:3
    for j=1:1:npts
        fpeBezierSample(j,i) = calcBezierYFcnXDerivative(lceN(j,1),...
                                              fiberForceLengthCurve,i-2);

        fpeTanhSample(j,i) = calcTanhSeriesDerivative(lceN(j,1),...
                                       forceLengthTanhCoeffs,i-2);
    end

    figure(figCurves);
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
            legend('Location','NorthEast');            
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

%%
% Make the force-velocity curve
%%

curvinessEccentricForceVelocity = 1.0;
flag_sharpEccentricTransition = 0;

forceVelocityMultiplierAtHalfMaximumFiberVelocity       = 0.2;
forceVelocityMultiplierAtLowEccentricFiberVelocity      = 1.4;
forceVelocityMultiplierAtMaximumEccentricFiberVelocity  = 1.6;

fiberForceVelocityCurve ...
  = createFiberForceVelocityCurve(...
      forceVelocityMultiplierAtHalfMaximumFiberVelocity,...
      forceVelocityMultiplierAtLowEccentricFiberVelocity,...
      forceVelocityMultiplierAtMaximumEccentricFiberVelocity,...
      curvinessEccentricForceVelocity,...
      flag_sharpEccentricTransition,...
      flag_enableNumericallyNonZeroGradients,...
      smallNumericallyNonZeroValue,...
      smallNumericallyNonZeroSlope,...
      'ForceVelocityCurve',...
      flag_usingOctave);



x0   = -1.05;
y0   = calcBezierYFcnXDerivative(x0,fiberForceVelocityCurve,0);
dydx0= 0;%calcBezierYFcnXDerivative(x0,fiberForceVelocityCurve,1);

x1   = 0;
y1   = calcBezierYFcnXDerivative(x1,fiberForceVelocityCurve,0);
dydx1= calcBezierYFcnXDerivative(x1,fiberForceVelocityCurve,1);


yNegInf = 0;
yInf    = inf;  

xScale      = 1.75;
xAtIntYZero = x0;
xPoint = x1;
yPoint = y1;

[A,B,C,D,E,F] = calcTanhSegmentCoefficientsUpd(...
                    x0,x1,dydx0,dydx1,...
                    yNegInf,yInf,...
                    xScale,xPoint, yPoint, xAtIntYZero);

%[A,B,C,D,E,F] = calcTanhSegmentCoefficients(x0,x1,dydx0,dydx1,...
%                                            yNegInf,yInf,...
%                                            xShift,xScale,xAtIntYZero);
forceVelocityTanhCoeffs = [A,B,C,D,E,F];  

npts = 100;
xP0 = -1.1;
xP1 = 1.3;
vceN = [(xP0):(xP1-xP0)/(npts-1):(xP1)]';

fvBezierSample = zeros(npts,3);
fvTanhSample   = zeros(npts,3);

for i=2:1:3
    for j=1:1:npts
        fvBezierSample(j,i) = calcBezierYFcnXDerivative(vceN(j,1),...
                                              fiberForceVelocityCurve,i-2);

        fvTanhSample(j,i) = calcTanhSeriesDerivative(vceN(j,1),...
                                       forceVelocityTanhCoeffs,i-2);
    end

    figure(figCurves);
    subplot('Position',reshape(subPlotPanel(3,i,:),1,4));

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
            plot([x0;x1],[y0;y1],'o','MarkerSize',5,...
                'Color',bezierColor,'MarkerFaceColor',[1,1,1]);
            hold on;
            title('Force-Velocity Curve Value');            
            xlabel('Norm. Velocity ($v/v^M_{max}$)')
            ylabel('Norm. Force ($f/f^M_o$)');            

            
        case 3
            plot([x0;x1],[dydx0;dydx1],'o','MarkerSize',5,...
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

rmpath(parametersDirectoryTreeMTParams);
rmpath(parametersDirectoryTreeExperiments);
rmpath(parametersDirectoryTreeModels);
rmpath(parametersDirectoryTreeCurves);

