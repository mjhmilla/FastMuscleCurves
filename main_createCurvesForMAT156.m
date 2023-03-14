

clc;
close all;
clear all;

opengl('save','software');

%%
% Architectural Parameters from the LS-DYNA files
%%
lceOpt      =  0.0428571;
fceOpt      = 22.469172;
penOpt      =  0.1221730;
penOptD     = 7.0000000;
lceOptAT    =  0.0425376;
fceOptAT    = 22.301690;
ltSlk       =  0.0304511;
et          =  0.045833;
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
load(fullfile('output','structs','defaultFelineSoleus.mat'));
load(fullfile('output','structs','oneFelineSoleus.mat'));
load(fullfile('output','structs','zeroFelineSoleus.mat'));

%Generate the normalized


%%
% Remove the directories
%%
rmpath(parametersDirectoryTreeMTParams);
rmpath(parametersDirectoryTreeExperiments);
rmpath(parametersDirectoryTreeModels);
rmpath(parametersDirectoryTreeCurves);
