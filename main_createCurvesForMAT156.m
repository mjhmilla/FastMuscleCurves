

clc;
close all;
clear all;

opengl('save','software');

%%
% Architectural Parameters from the LS-DYNA files
%%
% lceOpt      =  0.0428571;
% fceOpt      = 22.469172;
% penOpt      =  0.1221730;
% penOptD     = 7.0000000;
% lceOptAT    =  0.0425376;
% fceOptAT    = 22.301690;
% ltSlk       =  0.0304511;
% et          =  0.045833;
vceMax      = 4.5;

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

%%
% Add the directories needed to run this script
%%
parametersDirectoryTreeMTParams     = genpath('parameters');
parametersDirectoryTreeExperiments  = genpath('experiments');
parametersDirectoryTreeModels       = genpath('models');
parametersDirectoryTreeCurves       = genpath('curves');

addpath(parametersDirectoryTreeMTParams);
addpath(parametersDirectoryTreeExperiments);
addpath(parametersDirectoryTreeModels);
addpath(parametersDirectoryTreeCurves);

%%
% Load the curves that we've already generated for the cat soleus
%%

load('output/structs/defaultFelineSoleusQuadraticCurves.mat');

npts  = 100;
domain= [] ; %Take the default extended range

falValues = calcQuadraticBezierYFcnXCurveSampleVector(...
  felineSoleusNormMuscleQuadraticCurves.activeForceLengthCurve,...
  npts, domain);

fpeValues = calcQuadraticBezierYFcnXCurveSampleVector(...
  felineSoleusNormMuscleQuadraticCurves.fiberForceLengthCurve,...
  npts, domain);

fvValues = calcQuadraticBezierYFcnXCurveSampleVector(...
  felineSoleusNormMuscleQuadraticCurves.fiberForceVelocityCurve,...
  npts, domain);

fig=figure;
subplot(1,3,1);
    plot(falValues.x,falValues.y,'.','Color',[0,0,0]);
    hold on;
    xlabel('$$\tilde{\ell}^{M}$$');
    ylabel('$$\tilde{f}^{L}$$');
    box off;

subplot(1,3,2);
    plot(fpeValues.x,fpeValues.y,'.','Color',[0,0,0]);
    hold on;
    xlabel('$$\tilde{\ell}^{M}$$');
    ylabel('$$\tilde{f}^{PE}$$');
    box off;
    
subplot(1,3,3);
    plot(fvValues.x,fvValues.y,'.','Color',[0,0,0]);
    hold on;
    xlabel('$$\tilde{v}^{M}$$');
    ylabel('$$\tilde{f}^{V}$$');
    box off;


vN = [fvValues.x;1];
fvN = [fvValues.y;...
       calcQuadraticBezierYFcnXDerivative(1,...
       felineSoleusNormMuscleQuadraticCurves.fiberForceVelocityCurve,0)]; 

success = writeFortranVector(falValues.x, falValues.y, 10, ...
    'output/tables/FortranExport/MAT156Tables/defaultFelineSoleusQ_activeForceLengthCurve.f');
success = writeFortranVector(fpeValues.x, fpeValues.y, 11, ...
    'output/tables/FortranExport/MAT156Tables/defaultFelineSoleusQ_forceLengthCurve.f');
success = writeFortranVector(vN, fvN, 12, ...
    'output/tables/FortranExport/MAT156Tables/defaultFelineSoleusQ_forceVelocityCurve.f');


%%
% Remove the directories
%%
rmpath(parametersDirectoryTreeMTParams);
rmpath(parametersDirectoryTreeExperiments);
rmpath(parametersDirectoryTreeModels);
rmpath(parametersDirectoryTreeCurves);
