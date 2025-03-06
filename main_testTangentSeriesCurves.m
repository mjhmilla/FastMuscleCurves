clc;
close all;
clear all;

flag_plotActiveForceLengthCurves= 1;
flag_plotForceLengthCurves      = 0;
flag_plotForceVelocityCurves    = 0;
flag_plotTendonCurves           = 0;

flag_writePlotToFile=1;

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
bezierColor = [1,1,1].*0.85;
tanhColor   = [0,0,1];
tanColor    = [1,0,1];
expColor    = [0.5,0,0];

numberOfVerticalPlotRows      = 4;
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
indexPlotRow=1;
flag_usingOctave=0;

plotSettings.indexPlotRow=0;
plotSettings.subPlotPanel=subPlotPanel;
plotSettings.flag_plotBezierCurves=flag_plotBezierCurves;
plotSettings.flag_plotTanhCurves=flag_plotTanhCurves;
plotSettings.flag_plotTanCurves=flag_plotTanCurves;
plotSettings.bezierColor=bezierColor;
plotSettings.bezierWidth=2;
plotSettings.tanhColor=tanhColor;
plotSettings.tanColor=tanColor;
plotSettings.tanhErrorColor=tanhColor.*0.5 + [1,1,1].*0.5;
plotSettings.expColor = expColor;
plotSettings.expErrorColor = expColor.*0.5 + [1,1,1].*0.5;

%%
% Reference Bezier Curves
%%
flag_enableNumericallyNonZeroGradients  = 0;
smallNumericallyNonZeroNumber           = sqrt(sqrt(eps));
smallNumericallyNonZeroSlope            = sqrt(eps);
npts                                    = 101;

%
% Force velocity curve
%
curvinessEccentricForceVelocity = 1.0;
flag_sharpEccentricTransition = 0;

forceVelocityMultiplierAtHalfMaximumFiberVelocity       = 0.2;
forceVelocityMultiplierAtLowEccentricFiberVelocity      = 1.4;
forceVelocityMultiplierAtMaximumEccentricFiberVelocity  = 1.5;

fvCurveParams.fvAtHalfVceCMax = 0.2;
fvCurveParams.fvAtHalfVceELow = 1.4;
fvCurveParams.fvAtHalfVceEMax = 1.5;

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

vce0 = fiberForceVelocityCurve.xEnd(1,1);
vce1 = fiberForceVelocityCurve.xEnd(1,2);


% vceA = vce0;
% vceB = vceA + 0.1;
% fvCFast = [vceA:((vceB-vceA)/(9)):vceB];
% fvSlow = [-0.1:((0.2)/(9)):0.1];
% fvEFast = [0.9:(0.1/(9)):1.0];
%fvDomainTest = [fvCFast';fvSlow';fvEFast'];

fvDomainTest = [vce0:((vce1-vce0)/(npts-1)):vce1]';
%fvDomainTest = [fvDomainTest;vce1*1.1];

%
%Force-length curve
%

normLengthZero = 1+0; 
normLengthToe  = 1+0.6;
fToe  = 1;
yZero = 0;
kZero = 0;
if(flag_enableNumericallyNonZeroGradients==1)
    yZero   = smallNumericallyNonZeroNumber;
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

lce0 = fiberForceLengthCurve.xEnd(1,1);
lce1 = fiberForceLengthCurve.xEnd(1,2) ... 
     + 0.1*diff(fiberForceLengthCurve.xEnd);

fpeDomainTest = [lce0:((lce1-lce0)/(npts-1)):lce1]';
fpeDomainTest = [fpeDomainTest;(lce1+(lce1-lce0)*0.1)];

fpeCurveParams.lceZeroN      = normLengthZero;
fpeCurveParams.lceToeN       = normLengthToe;
fpeCurveParams.fToeN         = fToe;
fpeCurveParams.kToeN         = kToe;

%
%Active force length curve
%
assert(abs(fToe-1)<sqrt(eps));
normFiberLengthAtOneNormPassiveForce = normLengthToe;
normPevkToActinAttachmentPoint       = 0.5;
normMaxActiveTitinToActinDamping     = 65;
ecmForceFraction                     = 0.56;

[sarcomereProperties] =...
  getMammalianSkeletalMuscleNormalizedSarcomereProperties(...
    'human',...
    normFiberLengthAtOneNormPassiveForce,...
    normPevkToActinAttachmentPoint,...
    normMaxActiveTitinToActinDamping,...
    ecmForceFraction,...
    []);

curvinessActiveForceLength                       = 1.0;
shiftLengthActiveForceLengthCurveDescendingCurve = 0;
flag_useOctave = 0;

flag_compensateForCrossbridgeStiffness = 0;
[activeForceLengthCurve, ...
  activeForceLengthCurveAnnotationPoints] ...
    = createFiberActiveForceLengthCurve(...
          sarcomereProperties.normMyosinHalfLength*2,...
          sarcomereProperties.normMyosinBareHalfLength*2,...
          sarcomereProperties.normActinLength,...
          sarcomereProperties.normZLineLength,...
          sarcomereProperties.normSarcomereLengthZeroForce,...
          sarcomereProperties.normCrossBridgeStiffness,...          
          curvinessActiveForceLength, ...           
          shiftLengthActiveForceLengthCurveDescendingCurve,...
          flag_compensateForCrossbridgeStiffness,...          
          flag_enableNumericallyNonZeroGradients,...
          smallNumericallyNonZeroValue,...
          smallNumericallyNonZeroSlope,...
          'activeForceLength',...
          flag_useOctave);  

falCurveParams.scale = 1;
falCurveParams.x = activeForceLengthCurveAnnotationPoints.x(:,1);
falCurveParams.y = activeForceLengthCurveAnnotationPoints.y(:,1);

lceNMin = activeForceLengthCurve.xEnd(1,1) - 0.1;
lceNMax = activeForceLengthCurve.xEnd(1,2) + 0.1;

falDomainTest = [lceNMin:((lceNMax-lceNMin)/(npts-1)):lceNMax]';

%
% Tendon-force-length curve
%

ltsN            = 1.0;
eIso            = 0.049;
kIso            = 1.375/eIso;
fToe            = 2./3.;
curvinessTendon = 0.5;
computeIntegral = 1;
minimumSlope    = 0;
if(flag_enableNumericallyNonZeroGradients==1)
  minimumSlope = smallNumericallyNonZeroNumber/10.;
end

tendonForceLengthCurve = ...
  createTendonForceLengthCurve2021( ltsN, eIso, kIso, ...
                                    fToe, curvinessTendon, ...
                                    computeIntegral, ...
                                    flag_enableNumericallyNonZeroGradients,...
                                    smallNumericallyNonZeroNumber,...
                                    smallNumericallyNonZeroSlope,...
                                    'TendonForceLengthCurve',...
                                    flag_usingOctave);
ltN0 = ltsN;
ltN1 = ltsN+eIso*1.1;
npts = 100;
ftDomainTest = [ltN0:((ltN1-ltN0)/(npts-1)):ltN1]';

tendonCurveParams.ltZeroN = tendonForceLengthCurve.xEnd(1,1);
tendonCurveParams.ltToeN  = tendonForceLengthCurve.xEnd(1,2);
tendonCurveParams.fToeN   = tendonForceLengthCurve.yEnd(1,2);
tendonCurveParams.kToeN   = tendonForceLengthCurve.dydxEnd(1,2);
tendonCurveParams.ltIsoN  = eIso+1;

tendonCurveParamsDeGroote.eIso = eIso;

%%
% Tanh spline curve
%%
flag_useTanhSpline=1;

if(flag_useTanhSpline==1)



    xk = [activeForceLengthCurve.xpts(1,1:(end))';...
          activeForceLengthCurve.xpts(end,end)];
    yk      = zeros(size(xk));
    dydxk   = zeros(size(xk));
    for i=1:1:length(xk)
        yk(i,1)=calcBezierYFcnXDerivative(xk(i,1),activeForceLengthCurve,0);
        dydxk(i,1)=calcBezierYFcnXDerivative(xk(i,1),activeForceLengthCurve,1);
    end
    dydxk(1,1)=activeForceLengthCurve.dydxEnd(1,1);
    dydxk(end,1)=activeForceLengthCurve.dydxEnd(1,2);
    yLim = [0,0];
    
    nSample=100;
    x0 = activeForceLengthCurve.xEnd(1,1);
    x1 = activeForceLengthCurve.xEnd(1,2);
    xSample = [x0:((x1-x0)/(nSample-1)):x1]';
    xSample = sort([xSample;1]);
    ySample=zeros(size(xSample));
    for i=1:1:length(xSample)
        ySample(i,1) = calcBezierYFcnXDerivative(xSample(i,1),activeForceLengthCurve,0);
    end
    xAtIntYZero = 0;
    smoothness=1;
    flag_polishKnotPoints=1;
    indexOfValueKnotsToPolish = [1,4,6];
    indexOfDerivativeKnotsToPolish = [1,4,6];

    falTanhSplineCoeffs = ...
        fitTanhSplineCoefficients(xk,yk,dydxk,smoothness,...
            yLim,xAtIntYZero,flag_polishKnotPoints,...
            indexOfValueKnotsToPolish,indexOfDerivativeKnotsToPolish);
    
    falTanhSample = zeros(size(ySample));
    for i=1:1:length(xSample)
        falTanhSample(i,1)= ...
            calcTanhSeriesDerivative(xSample(i,1),falTanhSplineCoeffs,0);
    end

    flag_testPartialDerivatives=0;
    if(flag_testPartialDerivatives==1)
        fprintf(['Checking calcTanhSeriesParameterDerivative',...
                 ' against numerical derivatives.']);
        h=sqrt(eps);
        for derOrder=0:1
            fprintf('Checking partial derivatives of function derivative %d\n',derOrder);
            for i=1:1:6
                x=0.85;
                coeffs = falTanhSplineCoeffs(2,:);
                coeffsL = coeffs;
                coeffsL(1,i)=coeffsL(1,i)-h;
                coeffsR = coeffs;
                coeffsR(1,i)=coeffsR(1,i)+h;
                
                dP = calcTanhSeriesParameterDerivative(...
                        x,coeffs,derOrder,i);
                dPL = calcTanhSeriesDerivative(...
                        x,coeffsL,derOrder);
                dPR = calcTanhSeriesDerivative(...
                        x,coeffsR,derOrder);
                dPNum = (dPR-dPL)/(2*h);
    
                dPErr = abs(dP-dPNum);
                fprintf('%1d\t%1.3e\n',i,dPErr);
            end
        end
    end
end
%%
% Tanh curves
%%
if(flag_plotActiveForceLengthCurves==1)
    plotSettings.indexPlotRow=1;

    figCurves = addTanhActiveForceLengthCurveComparison(...
                    figCurves,...
                    falCurveParams,...
                    activeForceLengthCurve, ...
                    falDomainTest,...
                    plotSettings);

    if(flag_useTanhSpline==1)
        subplot('Position',...
            reshape(plotSettings.subPlotPanel(...
                plotSettings.indexPlotRow,2,:),1,4));
    
        plot(xSample,falTanhSample,'-m');
        hold on;
        plot(xk,yk,'xm');
        hold on;
    end

    plotSettings.flag_plotBezierCurves=0;
    figCurves = addDeGrooteFregly2016ActiveForceLengthCurveComparison(...
                    figCurves,...
                    falCurveParams,...
                    activeForceLengthCurve, ...
                    falDomainTest,...
                    plotSettings); 
    plotSettings.flag_plotBezierCurves=1;    
end

if(flag_plotForceVelocityCurves==1)
    plotSettings.indexPlotRow=2;

    figCurves = addTanhForceVelocityCurveComparison(...
                    figCurves,...
                    fiberForceVelocityCurve, ...
                    fvDomainTest,...
                    plotSettings);

    plotSettings.flag_plotBezierCurves=0;
    figCurves = addDeGrooteFregly2016ForceVelocityCurveComparison(...
                    figCurves,...
                    fvCurveParams,...
                    fiberForceVelocityCurve, ...
                    fvDomainTest,...
                    plotSettings); 
    plotSettings.flag_plotBezierCurves=1;
end

if(flag_plotForceLengthCurves==1)
    plotSettings.indexPlotRow=3;
    figCurves = addTanhForceLengthCurveComparison(...
                    figCurves,...
                    fpeCurveParams,...
                    fiberForceLengthCurve,...
                    fpeDomainTest,...
                    plotSettings);

    plotSettings.flag_plotBezierCurves=0;
    figCurves = addDeGrooteFregly2016PassiveForceLengthCurveComparison(...
                    figCurves,...
                    fpeCurveParams,...
                    fiberForceLengthCurve, ...
                    fpeDomainTest,...
                    plotSettings);
    plotSettings.flag_plotBezierCurves=1;
end

if(flag_plotTendonCurves==1)
    plotSettings.indexPlotRow=4;
    figCurves = addTanhTendonCurveComparison(...
                    figCurves,...
                    tendonCurveParams,...
                    tendonForceLengthCurve,...
                    ftDomainTest,...
                    plotSettings);
    
    plotSettings.flag_plotBezierCurves=0;
    figCurves = addDeGrooteFregly2016TendonCurveComparison(...
                    figCurves,...
                    tendonCurveParamsDeGroote,...
                    tendonForceLengthCurve,...
                    ftDomainTest,...
                    plotSettings);
    plotSettings.flag_plotBezierCurves=1;
end

set(figCurves,'Units','centimeters',...
    'PaperUnits','centimeters',...
    'PaperSize',[pageWidth pageHeight],...
    'PaperPositionMode','manual',...
    'PaperPosition',[0 0 pageWidth pageHeight]);     
set(figCurves,'renderer','painters');     
set(gcf,'InvertHardCopy','off')

if(flag_writePlotToFile==1)
    print('-dpdf', [pubOutputFolder,'fig_HyperbolicMuscleCurves.pdf']);
end

% rmpath(parametersDirectoryTreeMTParams);
% rmpath(parametersDirectoryTreeExperiments);
% rmpath(parametersDirectoryTreeModels);
% rmpath(parametersDirectoryTreeCurves);

