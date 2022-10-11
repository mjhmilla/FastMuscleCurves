
clc;
close all;
clear all;

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

%Number of subdivisions of the quintic Bezier curves to use when
%constructing the quadratic Bezier curves
numberOfQuadraticSubdivisions= 2;
numberOfCubicSubdivisions    = 1;
cubicZeroSecondDerivative    = 1;

fitCrossBridgeStiffnessDampingToKirch199490Hz=0;
flag_useFixedLambdaECM    = 0;

flag_plotEveryCurve       = 1;
flag_plotForceLengthDetail= 1;


pubOutputFolder                         = 'output/plots/MuscleCurves/';
postprocessingDirectoryTree             = genpath('postprocessing');
addpath(postprocessingDirectoryTree   );
addpath('colornames/');

[names,rgb] = colornames('SVG','FireBrick');

flag_embedECMandTitinFractionIntoCurves = 0;

flag_useOctave = 0;
flag_enableNumericallyNonZeroGradients    = 1;

flag_usingOctave    = 0;
plotWidth           = 6;
plotHeight          = 6;
pageWidth           = 21;
pageHeight          = 29.7;
plotHorizMarginCm   = 2;
plotVertMarginCm    = 2;
numberOfVerticalPlotRows      = 1;
numberOfHorizontalPlotColumns = 2;

plotConfigGeneric;

% 1. Active force length curve vs. data
% Solution: There were some initial descrepencies between the experimental force
%length data and a theoretical curve. These errors almost completely go
%away if it is assumed that the experimental recordings are of total 
%path length, rather than fiber length. In this case, when the elastiticy
%of the tendon is taken into account the theoretical active-force-length 
%curve and the transformed data nicely align.

%Failed attempt:
%This creates a cat soleus with an optimal fiber length of 58 mm: this
%is simply way too big to be realistic (given the data I'm working from)
flag_solveForOptimalFiberLengthOfBestFit  = 0; 

%Failed attempt:
shiftLengthActiveForceLengthCurveDescendingCurve = 0.;%...
%  (1/3)*( (1.154-1.087) + (1.23-1.162) + (1.077-1.039) );

smallNumericallyNonZeroNumber           = sqrt(sqrt(eps));

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
% Cat Soleus
%%
scaleOptimalFiberLengthCatSoleus        = 1.0; 
scaleMaximumIsometricTensionCatSoleus   = 1;



%%
% Human Soleus Model parameters
%%
scaleOptimalFiberLengthHumanSoleus      = 1.0; 
scaleMaximumIsometricTensionHumanSoleus = 1;




normMaxActiveTitinToActinDamping = 65;

normPevkToActinAttachmentPointZero = 0;
normPevkToActinAttachmentPointOne  = 1;
normPevkToActinAttachmentPointMid = 0.5;%^0.8175;
normPevkToActinAttachmentPointDefault=0.5;

normFiberLengthAtOneNormPassiveForceDefault = 1.367732948060934e+00;
% This the normalized fiber length at which the passive-force-length curve 
% used in this work develops 1 maximum isometric force passively when fit to 
% passive-force length data from Fig. 7 of Herzog & Leonard. This is used as 
% a default value for the other models used in this work for which the 
% passive force length properties are unknown.
%
% Herzog W, Leonard TR. Force enhancement following stretching of skeletal 
% muscle: a new mechanism. Journal of Experimental Biology. 
% 2002 May 1;205(9):1275-83.


ecmForceFractionDefault = 0.56;
% This is the average contribution of the ECM to the passive force length
% curve. Prado et al. reports the average contribution of titin to the 
% passive stiffness of 5 muscles in a rabbit. For a default value of the
% ECM's contribution we use 1-mean(titinContribution) = 0.56:
%
%Page 472 column 2 last paragraph of Prado et al. reports:
%
%     These results show that titin’s relative contribution to
%     total passive stiffness is much higher in some muscles,
%     like psoas (57%) and diaphragm (56%), than in oth-
%     ers, like soleus (24%), EDL (42%), and gastrocne-
%     mius (41%)
%
% Prado LG, Makarenko I, Andresen C, Krüger M, Opitz CA, Linke WA. Isoform 
% diversity of giant proteins in relation to passive and active contractile 
% properties of rabbit skeletal muscles. The Journal of general physiology. 
% 2005 Nov;126(5):461-80.


ecmForceFractionHumanSoleus     = ecmForceFractionDefault;
ecmForceFractionFelineSoleus    = ecmForceFractionDefault;

wlcTitinModel = 1;
% This titin model type will extend the force length curve of each titin
% segment to fit the WLC model up to a high force without going to a 
% singularity.

linearTitinModel=0;
% This titin model type will linear extrapolate the force length curve
% from the length at which the muscle develops a passive force equivalent
% to one maximum sometric force.

useCalibratedCurves     = 1;
% Here calibrated curves refers to an active-force-length curve and a 
% force velocity curve that have been adjusted so that opus 31 (the
% proposed model) can reproduce the desired active-force-length and force
% velocity curves. This is necessary because the deformation of the 
% viscoelastic cross-bridge element is not taken into consideration 
% in the formulation of these curves, but does affect the output.


useTwoSidedTitinCurves  = 0;
% useTwoSidedTitinCurves = 0? 
% The force length curve of the prox. IG and PEVK segments has a slack 
% length: below the slack length the force goes to zero just like a tendon.
% 
% useTwoSidedTitinCurves = 1?
% The elastic titin curves have a shape, broadly speaking, like a tan 
% function. In the case, that the titin-actin bond is formed at a long 
% sarcomere length and then the sarcomere is rapidly shortened, the distal
% section of titin may become more proximal to the Z-line than the 
% titin-actin bond (a figure is really needed to properly explain this):
% in these circumstances the PEVK segment would have a negative length
% and would generate an elastic forces that impedes further shortening.
% If this came to pass, titin would reduce the tension a muscle can generate
% during shortening. I found this not to be the case during the ramp 
% shortening simulations used in this work and so by default I set 
% useTwoSidedTitinCurves=0. It migth be the case that a shortening ramp 
% beginning at a longer length would see some additional force reduction
% due to this effect.

titinModel = linearTitinModel; %wlcTitinModel;

[ defaultFelineSoleus,...
  activeForceLengthCurveAnnotationPoints,...
  felineSoleusActiveForceLengthDataDefault,...
  felineSoleusPassiveForceLengthDataDefault,...
  felineSoleusPassiveForceLengthCurveSettings ] = ...
        createFelineSoleusModel(...
                normPevkToActinAttachmentPointDefault,...
                normMaxActiveTitinToActinDamping,...
                normFiberLengthAtOneNormPassiveForceDefault,... 
                ecmForceFractionFelineSoleus,...
                titinModel,...
                useCalibratedCurves,...
                useTwoSidedTitinCurves,...
                smallNumericallyNonZeroNumber,...
                flag_enableNumericallyNonZeroGradients,...
                scaleOptimalFiberLengthCatSoleus,... 
                scaleMaximumIsometricTensionCatSoleus,...
                flag_embedECMandTitinFractionIntoCurves,...
                flag_useOctave);

save('output/structs/defaultFelineSoleus.mat',...
     'defaultFelineSoleus');  


zeroHumanSoleus = createHumanSoleusModel(...
                        normPevkToActinAttachmentPointZero,...
                        normMaxActiveTitinToActinDamping,...                        
                        normFiberLengthAtOneNormPassiveForceDefault,... 
                        ecmForceFractionHumanSoleus,...
                        titinModel,...
                        useCalibratedCurves,...
                        useTwoSidedTitinCurves,...
                        smallNumericallyNonZeroNumber,...
                        flag_enableNumericallyNonZeroGradients,...
                        scaleOptimalFiberLengthHumanSoleus,...
                        scaleMaximumIsometricTensionHumanSoleus,...
                        felineSoleusPassiveForceLengthCurveSettings,...
                        flag_embedECMandTitinFractionIntoCurves,...
                        flag_useOctave);

oneHumanSoleus = createHumanSoleusModel(...
                        normPevkToActinAttachmentPointOne,...
                        normMaxActiveTitinToActinDamping,...                        
                        normFiberLengthAtOneNormPassiveForceDefault,... 
                        ecmForceFractionHumanSoleus,...
                        titinModel,...
                        useCalibratedCurves,...
                        useTwoSidedTitinCurves,...
                        smallNumericallyNonZeroNumber,...
                        flag_enableNumericallyNonZeroGradients,...
                        scaleOptimalFiberLengthHumanSoleus,...
                        scaleMaximumIsometricTensionHumanSoleus,...
                        felineSoleusPassiveForceLengthCurveSettings,...
                        flag_embedECMandTitinFractionIntoCurves,...
                        flag_useOctave);

midHumanSoleus = createHumanSoleusModel(...
                        normPevkToActinAttachmentPointMid,...
                        normMaxActiveTitinToActinDamping,...                        
                        normFiberLengthAtOneNormPassiveForceDefault,... 
                        ecmForceFractionHumanSoleus,...
                        titinModel,...
                        useCalibratedCurves,...
                        useTwoSidedTitinCurves,...
                        smallNumericallyNonZeroNumber,...
                        flag_enableNumericallyNonZeroGradients,...
                        scaleOptimalFiberLengthHumanSoleus,...
                        scaleMaximumIsometricTensionHumanSoleus,...
                        felineSoleusPassiveForceLengthCurveSettings,...
                        flag_embedECMandTitinFractionIntoCurves,...
                        flag_useOctave);



figInterpTitin = figure;

n01 = [0:0.01:1]';

structNamesPlot = {'zeroHumanSoleus',...
                   'midHumanSoleus',...
                   'oneHumanSoleus'};

structColors = [1,0,0;...
                1,0,1;...
                0,0,1];

curveNamesPlot  = {'forceLengthProximalTitinCurve',...
                   'forceLengthDistalTitinCurve'};

xSamples = zeros(length(n01),1);
ySamples = zeros(length(n01),1);

%Plot the interpolated curves

maxBlend = 11;

for idxCurve = 1:1:length(curveNamesPlot)
    
    for idxBlend = 1:1:maxBlend

        curveStruct = [];

        curveStructA = zeroHumanSoleus.curves.(curveNamesPlot{idxCurve});
        curveStructB = oneHumanSoleus.curves.(curveNamesPlot{idxCurve});

        A = 1 - ((idxBlend-1)/(maxBlend-1));
        B = 1 - A;

        curveStruct = curveStructA;

        curveStruct.name = sprintf('%d',round(A*100,2));
        curveStruct.xpts = A.*curveStructA.xpts + B.*curveStructB.xpts;
        curveStruct.ypts = A.*curveStructA.ypts + B.*curveStructB.ypts;
        curveStruct.xEnd = A.*curveStructA.xEnd + B.*curveStructB.xEnd;
        curveStruct.yEnd = A.*curveStructA.yEnd + B.*curveStructB.yEnd;
        curveStruct.dydxEnd = [0,0].*nan;
        curveStruct.d2ydx2End = [0,0].*nan;
        curveStruct.integral = [];

        if(idxCurve==1 && idxBlend==1)
            disp('First and second derivatives cannot be averaged at the ends');
            disp('... I will probably have to calculate and cache these values');
            disp('... and only update if the parameter is changed');
        end

        dydx0 = calcBezierYFcnXDerivative(curveStruct.xEnd(1,1), curveStruct, 1);
        dydx1 = calcBezierYFcnXDerivative(curveStruct.xEnd(1,2), curveStruct, 1);

        d2ydx20 = calcBezierYFcnXDerivative(curveStruct.xEnd(1,1), curveStruct, 2);
        d2ydx21 = calcBezierYFcnXDerivative(curveStruct.xEnd(1,2), curveStruct, 2);
    
        curveStruct.dydxEnd = [dydx0,dydx1];
        curveStruct.d2ydx2End = [d2ydx20,d2ydx21];
        

        xEnd = curveStruct.xEnd;
        xA = xEnd(1,1)-0.1*(diff(xEnd));
        xB = xEnd(1,2)+0.1*(diff(xEnd));        
    
        for idxPt=1:1:length(n01)            
    
            xS = xA + n01(idxPt,1)*(xB-xA);
            yS = calcBezierYFcnXDerivative(xS, curveStruct, 0);
    
            xSamples(idxPt,1) = xS;
            ySamples(idxPt,1) = yS;
                   
        end
    
        lineColor = structColors(1,:).*A + structColors(end,:).*B;

        subplot('Position', reshape(subPlotPanel(1,idxCurve,:),1,4));
    
        plot(xSamples(:,1),ySamples(:,1),...
            'Color',lineColor.*0.5 + [1,1,1].*0.5,...
            'LineWidth', 2);
        hold on;
        text(xSamples(end,1),ySamples(end,1)+0.05,curveStruct.name,...
             'Color',lineColor.*0.5 + [1,1,1].*0.5);
        hold on;
    end
    
    title(curveNamesPlot{idxCurve});
    xlabel('Norm. Length');
    ylabel('Norm. Force');
    
    box off;
    %legend;
    %legend boxoff;

end


%Plot the constructed curves
for idxStruct = 1:1:length(structNamesPlot)
    for idxCurve = 1:1:length(curveNamesPlot)
        
        curveStruct = [];
        seriesLabel = '';
        if(contains(structNamesPlot{idxStruct},'zeroHumanSoleus'))
            curveStruct = zeroHumanSoleus.curves.(curveNamesPlot{idxCurve});
            seriesLabel = '0';
        end
        if(contains(structNamesPlot{idxStruct},'midHumanSoleus'))
            curveStruct = midHumanSoleus.curves.(curveNamesPlot{idxCurve});
            seriesLabel = '0.5';
        end
        if(contains(structNamesPlot{idxStruct},'oneHumanSoleus'))
            curveStruct = oneHumanSoleus.curves.(curveNamesPlot{idxCurve});
            seriesLabel = '1';
        end

        idxData = (idxStruct-1)*length(curveNamesPlot) + idxCurve;

        xEnd = curveStruct.xEnd;
        xA = xEnd(1,1)-0.1*(diff(xEnd));
        xB = xEnd(1,2)+0.1*(diff(xEnd));        

        for idxPt=1:1:length(n01)            

            xS = xA + n01(idxPt,1)*(xB-xA);
            yS = calcBezierYFcnXDerivative(xS, curveStruct, 0);

            xSamples(idxPt,1) = xS;
            ySamples(idxPt,1) = yS;
                   
        end

        subplot('Position', reshape(subPlotPanel(1,idxCurve,:),1,4));

        plot(xSamples(:,1),ySamples(:,1),...
            'Color',[0,0,0],...
            'LineWidth', 1,...
            'DisplayName',structNamesPlot{idxStruct});
        hold on;

        text(xSamples(end,1),ySamples(end,1),seriesLabel);
        hold on;

        if(idxData==1 || ...
           idxData == length(curveNamesPlot)*length(structNamesPlot))
            
            title(curveNamesPlot{idxCurve});
            xlabel('Norm. Length');
            ylabel('Norm. Force');
            
            box off;
            %legend;
            %legend boxoff;
        end
    end

end


curveNames  = {'forceLengthProximalTitinCurve',...
               'forceLengthDistalTitinCurve'};

xError = zeros(size(zeroHumanSoleus.curves.(curveNames{1}).xpts));
yError = zeros(size(zeroHumanSoleus.curves.(curveNames{1}).ypts));

for c=1:1:length(curveNames)

    for i=1:1:size(zeroHumanSoleus.curves.(curveNames{c}).xpts,1)
        for j=1:1:size(zeroHumanSoleus.curves.(curveNames{c}).xpts,2)

            a = normPevkToActinAttachmentPointMid;
            b = (1-a);

            xInt =  (b*zeroHumanSoleus.curves.(curveNames{c}).xpts(i,j)) ...
                   +(a *oneHumanSoleus.curves.(curveNames{c}).xpts(i,j));
            xMid = midHumanSoleus.curves.(curveNames{c}).xpts(i,j);
            xErr = xInt-xMid;
            xError(i,j)=xErr;
    
            yInt =  (b*zeroHumanSoleus.curves.(curveNames{c}).ypts(i,j)) ...
                   +(a *oneHumanSoleus.curves.(curveNames{c}).ypts(i,j));
            yMid = midHumanSoleus.curves.(curveNames{c}).ypts(i,j);
            yErr  = yInt-yMid;
            yError(i,j)=yErr;    
    
        end
    end
    disp(curveNames{c});
    disp('x-error');
    disp(xError);
    disp('y-error');
    disp(yError);
    
end


