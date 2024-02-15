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
indexPlotRow=1;
flag_usingOctave=0;

plotSettings.indexPlotRow=0;
plotSettings.subPlotPanel=subPlotPanel;
plotSettings.flag_plotBezierCurves=flag_plotBezierCurves;
plotSettings.flag_plotTanhCurves=flag_plotTanhCurves;
plotSettings.flag_plotTanCurves=flag_plotTanCurves;
plotSettings.bezierColor=bezierColor;
plotSettings.tanhColor=tanhColor;
plotSettings.tanColor=tanColor;

plotSettings.indexPlotRow=1;
figCurves = addTanhForceVelocityCurveComparison(figCurves,plotSettings,...
                                flag_usingOctave);

plotSettings.indexPlotRow=2;
figCurves = addTanhForceLengthCurveComparison(figCurves,plotSettings,...
                                flag_usingOctave);

plotSettings.indexPlotRow=3;
figCurves = addTanhTendonCurveComparison(figCurves,plotSettings,...
                                flag_usingOctave);






rmpath(parametersDirectoryTreeMTParams);
rmpath(parametersDirectoryTreeExperiments);
rmpath(parametersDirectoryTreeModels);
rmpath(parametersDirectoryTreeCurves);

