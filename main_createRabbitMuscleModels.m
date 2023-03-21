clc;
close all;
clear all;

parametersDirectoryTreeCurves       = genpath('curves');
addpath(parametersDirectoryTreeCurves);
postprocessingDirectoryTree         = genpath('postprocessing');
addpath(postprocessingDirectoryTree);


%%
% The purpose of this file is to create rabbit TA and EDL musculoskeletal models
% to simulate the experiments of
%
% Hasselman CT, Best TM, Seaber AV, Garrett JR WE. A threshold and continuum of 
% injury during active stretch of rabbit skeletal muscle. The American Journal 
% of Sports Medicine. 1995 Jan;23(1):65-73.
%
% The most complete measurements of the properties of rabbit TA and EDL leg 
% muscles (that I'm aware of) comes from
%
% Siebert T, Leichsenring K, Rode C, Wick C, Stutzig N, Schubert H, Blickhan R, 
% Böl M. Three-dimensional muscle architecture and comprehensive dynamic 
% properties of rabbit gastrocnemius, plantaris and soleus: input for 
% simulation studies. PLoS one. 2015 Jun 26;10(6):e0130985.
%
% Lieber RL, Blevins FT. Skeletal muscle architecture of the rabbit hindlimb: 
% functional implications of muscle design. Journal of morphology. 
% 1989 Jan;199(1):93-101.
%
% Here I'm going to implement the active force-length, passive-force-length,
% force-velocity, and tendon-force-length curves described in Siebert et al.
% and then identify the corresponding parameters for the model that I've
% developed.
%%

fal = struct('l1N',0,'l2N',0,'l3N',0,'l4N',0,'fc',0,'Fim',0,'lccopt',0);
fv  = struct('vccmax',0,'curv',0);
ft  = struct('F1_Fim',0,'lsec1_lsec0',0,'ksh',0,'k',0,'lsec0',0);
fpe = struct('k1',0,'k2',0,'lpec0',0);

edl.fal = fal;
edl.fv  = fv;
edl.ft  = ft;
edl.fpe = fpe;

ta.fal  = fal;
ta.fv   = fv;
ta.ft   = ft;
ta.fpe  = fpe;

%From Siebert et al. Supp. material
ta.fal.l1N      = 0.40;
ta.fal.l2N      = 0.80;
ta.fal.l3N      = 1.16;
ta.fal.l4N      = 2.45;
ta.fal.fc       = 0.84;
ta.fal.Fim      = 21.9;
ta.fal.lccopt   = 36.7;

edl.fal.l1N     = 0.39;
edl.fal.l2N     = 0.82;
edl.fal.l3N     = 1.21;
edl.fal.l4N     = 2.45;
edl.fal.fc      = 0.91;
edl.fal.Fim     = 13.0;
edl.fal.lccopt  = 14.1;

ta.fv.vccmax    = 16.4;
ta.fv.curv      = 0.33;

edl.fv.vccmax   = 12.3;
edl.fv.curv     = 0.37;

%The eccentric side of the curve is not described in Siebert et al. 2015
%Here we assume that the eccentric side is given by
%1+normEccentricDamping*vceN
normEccentricDamping = 0.30;

ta.ft.F1_Fim        = 0.20;
ta.ft.lsec1_lsec0   = 0.027;
ta.ft.ksh           = 2.1;
ta.ft.k             = 9.0;
ta.ft.lsec0         = 56.2;

edl.ft.F1_Fim        = 0.41;
edl.ft.lsec1_lsec0   = 0.035;
edl.ft.ksh           = 2.0;
edl.ft.k             = 15.8;
edl.ft.lsec0         = 90.3;

ta.fpe.k1       =  0.003; 
ta.fpe.k2       =  0.49; 
ta.fpe.lpec0    = 35.7; 
ta.fpe.Fim      = ta.fal.Fim;

edl.fpe.k1      =  0.001; 
edl.fpe.k2      =  0.82; 
edl.fpe.lpec0   = 13.7; 
edl.fpe.Fim     = edl.fal.Fim;

%%
% Sample the curves
%%

lceNRangeA = [0.3:0.01:3.0]';
lceNRangeB = [0.9:0.01:1.5]';
vceNRange = [-1.1:0.01:1.1]';
ltNRange  = [0.975:0.001:1.07]';


taSamples.falN  = zeros(size(lceNRangeA));
taSamples.fpeN  = zeros(size(lceNRangeB));
taSamples.fvN   = zeros(size(vceNRange));
taSamples.ftN   = zeros(size(ltNRange));

edlSamples.falN  = zeros(size(lceNRangeA));
edlSamples.fpeN  = zeros(size(lceNRangeB));
edlSamples.fvN   = zeros(size(vceNRange));
edlSamples.ftN   = zeros(size(ltNRange));

%Force-length properties
for i=1:1:length(lceNRangeA)
    lengthInMM = lceNRangeA(i,1)*ta.fal.lccopt;
    taSamples.falN(i,1) =...
        calcSiebert2015ActiveForceLengthCurve(lengthInMM, ta.fal);
end
for i=1:1:length(lceNRangeB)
    lengthInMM = lceNRangeB(i,1)*ta.fpe.lpec0;    
    taSamples.fpeN(i,1) = ...
        calcSiebert2015ForceLengthCurve(lengthInMM, ta.fpe);

end

%Force-velocity properties
for i=1:1:length(vceNRange)
    velInLPS = vceNRange(i,1)*ta.fv.vccmax;
    taSamples.fvN(i,1) = ...
        calcSiebert2015ForceVelocityCurve( ...
            velInLPS, ta.fv, normEccentricDamping);
end

%Force-length tendon properties
for i=1:1:length(ltNRange)
    lengthInMM = ta.ft.lsec0*ltNRange(i,1);

    taSamples.ftN(i,1) = ...
        calcSiebert2015TendonForceLengthCurve(lengthInMM, ta.ft);

end

%%
% Plotting settings
%%

opengl('save','software');

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

plotWidth = 6;
plotHeight= 6;
plotHorizMarginCm = 1.5;
plotVertMarginCm  = 2.0;

numberOfHorizontalPlotColumnsGeneric = 4;
numberOfVerticalPlotRowsGeneric      = 2;
                              
[subPlotPanelGeneric, pageWidthGeneric,pageHeightGeneric]= ...
      plotConfigGeneric(  numberOfHorizontalPlotColumnsGeneric,...
                          numberOfVerticalPlotRowsGeneric,...
                          plotWidth,plotHeight,...
                          plotHorizMarginCm,plotVertMarginCm); 

edlColor = [0,0,1]*0.5;
taColor  = [0,0,1];

fig=figure;

subplot('Position',reshape(subPlotPanelGeneric(1,1,:),1,4));
    plot(lceNRangeA, taSamples.falN,'Color',taColor,'DisplayName','TA');
    hold on;    

    legend;
    legend boxoff;
    box off;
    xlabel('Norm. Length $$\ell/\ell^{M}_o$$');
    ylabel('Norm. Force $$f/f^{M}_o$$');
    title('Active Force-Length');

subplot('Position',reshape(subPlotPanelGeneric(1,2,:),1,4));
    plot(lceNRangeB-1, taSamples.fpeN,'Color',taColor);
    hold on;    
    box off;
    xlabel('Norm. Length $$\ell/\ell^{PE}_o$$');
    ylabel('Norm. Force $$f/f^{M}_o$$');
    title('Passive Force-Length');
    xlim([0,1]);
    ylim([0,0.2]);

subplot('Position',reshape(subPlotPanelGeneric(1,3,:),1,4));
    plot(vceNRange, taSamples.fvN,'Color',taColor);
    hold on;    
    box off;
    xlabel('Norm. Velocity $$v/v^{M}_{max}$$');
    ylabel('Norm. Force $$f/f^{M}_o$$');
    title('Force Velocity Relation');

subplot('Position',reshape(subPlotPanelGeneric(1,4,:),1,4));
    plot(ltNRange, taSamples.ftN,'Color',taColor);
    hold on;    
    box off;    
    xlabel('Norm. Length $$\ell/\ell^{T}_{s}$$');
    ylabel('Norm. Force $$f/f^{M}_o$$');
    title('Tendon Force-Length');





