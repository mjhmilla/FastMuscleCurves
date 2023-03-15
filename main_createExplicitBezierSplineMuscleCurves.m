

clc;
close all;
clear all;

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

%Number of subdivisions of the quintic Bezier curves to use when
%constructing the quadratic Bezier curves
numberOfQuadraticSubdivisions= 2;
numberOfCubicSubdivisions    = 1;
cubicZeroSecondDerivative    = 1;

fitCrossBridgeStiffnessDampingToKirch199490Hz=0;
flag_useFixedLambdaECM    = 0;

flag_plotEveryCurve       = 1;
flag_plotForceLengthDetail= 1;

flag_plotQuadraticFelineSoleus  =1;
flag_plotQuadraticHumanSoleus   =1;

pubOutputFolder                         = 'output/plots/MuscleCurves/';
postprocessingDirectoryTree             = genpath('postprocessing');
addpath(postprocessingDirectoryTree   );
addpath('colornames/');

[names,rgb] = colornames('SVG','FireBrick');

flag_embedECMandTitinFractionIntoCurves = 0;

flag_useOctave = 0;
flag_enableNumericallyNonZeroGradients    = 1;

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
scaleOptimalFiberLengthCatSoleus        = 1; 
scaleMaximumIsometricTensionCatSoleus   = 1;
shiftExperimentalForceLengthDataCatSoleus = 0;

%%
% Human Soleus Model parameters
%%
scaleOptimalFiberLengthHumanSoleus      = 1.0; 
scaleMaximumIsometricTensionHumanSoleus = 1;



normMaxActiveTitinToActinDamping = 65;


normPevkToActinAttachmentPointDefault = 0.5;
normPevkToActinAttachmentPointZero = 0;
normPevkToActinAttachmentPointOne  = 1;
%The default value for the point of attachment between the PEVK segment
%and actin. This point of attachment is expressed as a normalized length
%where 0 corresponds to the start of the PEVK segment (at the prox Ig/PEVK
% boundary), 0.5 would be the middle of the PEVK segment, and 1.0 would 
% be the distal end of the PEVK segment (at the PEVK/distal Ig border).


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


fprintf(['\n\ndefaultFelineSoleus\n\n']);


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
                linearTitinModel,...
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

[ zeroFelineSoleus,...
  zeroActiveForceLengthCurveAnnotationPoints,...
  zeroFelineSoleusActiveForceLengthDataDefault,...
  zeroFelineSoleusPassiveForceLengthDataDefault,...
  zeroFelineSoleusPassiveForceLengthCurveSettings ] = ...
        createFelineSoleusModel(...
                normPevkToActinAttachmentPointZero,...
                normMaxActiveTitinToActinDamping,...
                normFiberLengthAtOneNormPassiveForceDefault,... 
                ecmForceFractionFelineSoleus,...
                linearTitinModel,...
                useCalibratedCurves,...
                useTwoSidedTitinCurves,...
                smallNumericallyNonZeroNumber,...
                flag_enableNumericallyNonZeroGradients,...
                scaleOptimalFiberLengthCatSoleus,... 
                scaleMaximumIsometricTensionCatSoleus,...
                flag_embedECMandTitinFractionIntoCurves,...
                flag_useOctave);

save('output/structs/zeroFelineSoleus.mat',...
     'zeroFelineSoleus');  


[ oneFelineSoleus,...
  oneActiveForceLengthCurveAnnotationPoints,...
  oneFelineSoleusActiveForceLengthDataDefault,...
  oneFelineSoleusPassiveForceLengthDataDefault,...
  oneFelineSoleusPassiveForceLengthCurveSettings ] = ...
        createFelineSoleusModel(...
                normPevkToActinAttachmentPointOne,...
                normMaxActiveTitinToActinDamping,...
                normFiberLengthAtOneNormPassiveForceDefault,... 
                ecmForceFractionFelineSoleus,...
                linearTitinModel,...
                useCalibratedCurves,...
                useTwoSidedTitinCurves,...
                smallNumericallyNonZeroNumber,...
                flag_enableNumericallyNonZeroGradients,...
                scaleOptimalFiberLengthCatSoleus,... 
                scaleMaximumIsometricTensionCatSoleus,...
                flag_embedECMandTitinFractionIntoCurves,...
                flag_useOctave);

save('output/structs/oneFelineSoleus.mat',...
     'oneFelineSoleus');


fprintf(['\n\ndefaultHumanSoleus\n\n']);
defaultHumanSoleus = createHumanSoleusModel(...
                        normPevkToActinAttachmentPointDefault,...
                        normMaxActiveTitinToActinDamping,...                        
                        normFiberLengthAtOneNormPassiveForceDefault,... 
                        ecmForceFractionHumanSoleus,...
                        linearTitinModel,...
                        useCalibratedCurves,...
                        useTwoSidedTitinCurves,...
                        smallNumericallyNonZeroNumber,...
                        flag_enableNumericallyNonZeroGradients,...
                        scaleOptimalFiberLengthHumanSoleus,...
                        scaleMaximumIsometricTensionHumanSoleus,...
                        felineSoleusPassiveForceLengthCurveSettings,...
                        flag_embedECMandTitinFractionIntoCurves,...
                        flag_useOctave);


save('output/structs/defaultHumanSoleus.mat',...
     'defaultHumanSoleus');                      




fprintf(['\n\nzeroHumanSoleus\n\n']);
zeroHumanSoleus = createHumanSoleusModel(...
                        normPevkToActinAttachmentPointZero,...
                        normMaxActiveTitinToActinDamping,...                        
                        normFiberLengthAtOneNormPassiveForceDefault,... 
                        ecmForceFractionHumanSoleus,...
                        linearTitinModel,...
                        useCalibratedCurves,...
                        useTwoSidedTitinCurves,...
                        smallNumericallyNonZeroNumber,...
                        flag_enableNumericallyNonZeroGradients,...
                        scaleOptimalFiberLengthHumanSoleus,...
                        scaleMaximumIsometricTensionHumanSoleus,...
                        felineSoleusPassiveForceLengthCurveSettings,...
                        flag_embedECMandTitinFractionIntoCurves,...
                        flag_useOctave);


fprintf(['\n\noneHumanSoleus\n\n']);
oneHumanSoleus = createHumanSoleusModel(...
                        normPevkToActinAttachmentPointOne,...
                        normMaxActiveTitinToActinDamping,...                        
                        normFiberLengthAtOneNormPassiveForceDefault,... 
                        ecmForceFractionHumanSoleus,...
                        linearTitinModel,...
                        useCalibratedCurves,...
                        useTwoSidedTitinCurves,...
                        smallNumericallyNonZeroNumber,...
                        flag_enableNumericallyNonZeroGradients,...
                        scaleOptimalFiberLengthHumanSoleus,...
                        scaleMaximumIsometricTensionHumanSoleus,...
                        felineSoleusPassiveForceLengthCurveSettings,...
                        flag_embedECMandTitinFractionIntoCurves,...
                        flag_useOctave);




%%
%Plot feline curves + experimental data from Herzog & Leonard 2002
%%

structOfFigures = [];

if(flag_plotEveryCurve==1)
    figH = plotStructOfBezierSplines( defaultHumanSoleus.curves,...
                                      {'Inverse','use'});                          
    
    %%
    % Note the average offset between the active-force-length curve and
    % the transformed data
    %%
    
    xExp = felineSoleusActiveForceLengthDataDefault(2:end,1);
    yExp = felineSoleusActiveForceLengthDataDefault(2:end,2);
    xCurve = zeros(size(xExp));
    
    for i=1:1:length(xExp)
    xCurve(i,1) = calcBezierFcnXGivenY(yExp(i,1), ...
      defaultFelineSoleus.curves.activeForceLengthCurve,... 
      xExp(i,1));
    end                                    
    dx = mean(xCurve-xExp);
    %felineSoleusActiveForceLengthDataDefault(:,1)=...
    %  felineSoleusActiveForceLengthDataDefault(:,1)+dx;
    
    disp('Normalized length offset');
    fprintf('%1.6f lce/lopt \n',felineSoleusActiveForceLengthDataDefault(1,1));
    disp('Average error on the descending limb');
    fprintf('%1.6f lce/lopt \n',dx);
    fprintf('%1.6f mm \n',dx*(defaultFelineSoleus.musculotendon.optimalFiberLength*1000));
    
    lceNStart = felineSoleusActiveForceLengthDataDefault(1,1);
    save('output/structs/normalizedFiberLengthStartHerzogLeonard2002.mat',...
         'lceNStart');
    
    % dl = felineSoleusActiveForceLengthDataDefault(2:end,1)-1;
    % A  = [dl ones(size(dl))];
    % b  = felineSoleusActiveForceLengthDataDefault(2:end,2);
    % 
    % x     = (A'*A)\(A'*b);
    % y0    = x(2,1);
    % dydx0 = x(1,1);
    % 
    % felineSoleusActiveForceLengthLineBestFit = zeros(length(dl),1);
    % 
    % felineSoleusActiveForceLengthLineBestFit(:,1) = ...
    %   felineSoleusActiveForceLengthDataDefault(2:end,1);
    % dl = felineSoleusActiveForceLengthLineBestFit(:,1)-1;
    % 
    % felineSoleusActiveForceLengthLineBestFit(:,2) = dydx0.*dl + y0;
    % 
    % disp('Active force length line of best fit y=(dydx)*(x-1) + y0');
    % fprintf('dydx: %1.3f\n',dydx0);
    % fprintf('  y0: %1.3f\n',y0);
    
    
    
    %%
    % Plot the derivative of the tendon force length curve on top of the
    % stiffness curve
    %% 
    figure(figH.tendonStiffnessCurve);
    
        curveSample = calcBezierYFcnXCurveSampleVector(...
                        defaultFelineSoleus.curves.('tendonForceLengthCurve'), 200,[]);
    
        xmin = min(curveSample.x);
        xmax = max(curveSample.x);
        ymin = min(curveSample.y);
        ymax = max(curveSample.y);
    
    subplot(2,2,1);
      plot(curveSample.x, curveSample.dydx,...
        '--','Color',[1,1,1].*0.5,'LineWidth',2);
      hold on;
    
    %%
    %Plot experimental data over top of the curves where it is available.
    %%
    figure(figH.activeForceLengthCurve);
      subplot(2,2,1);  
      plot(  felineSoleusActiveForceLengthDataDefault(:,1)+shiftExperimentalForceLengthDataCatSoleus,...
           felineSoleusActiveForceLengthDataDefault(:,2),'xb');
      hold on;          
    
    figure(figH.fiberForceLengthCurve);
      subplot(2,2,1);
      plot(   felineSoleusPassiveForceLengthDataDefault(:,1)+shiftExperimentalForceLengthDataCatSoleus,...
              felineSoleusPassiveForceLengthDataDefault(:,2),'xb');
      hold on;          
  
end

%%
% Plot the passive force length curve of the 2 segment titin model and
% compare it to the passive force length curve
%%
if(flag_plotForceLengthDetail==1)

    fpeDomain = defaultHumanSoleus.curves.('fiberForceLengthCurve').xEnd ...
                 +[-0.01,0.01]+[0,1];  
    n = 100;
    samples = [min(fpeDomain):((max(fpeDomain)-min(fpeDomain))/(n-1)):max(fpeDomain)]';
    
    [ curveSampleECMHalfHuman,...
      curveSampleTitinHuman,...
      curveSampleTitinActiveHuman,...
      curveSampleIgpHuman,...
      curveSamplePevkHuman,...
      curveSampleIgdHuman,...
      curveSampleProximalTitinHuman,...
      curveSampleDistalTitinHuman] = ...
      sampleTitinCurves(samples.*0.5,...
                        defaultHumanSoleus.curves.forceLengthECMHalfCurve,...
                        defaultHumanSoleus.curves.forceLengthProximalTitinCurve,...
                        defaultHumanSoleus.curves.forceLengthProximalTitinInverseCurve,...
                        defaultHumanSoleus.curves.forceLengthDistalTitinCurve,...
                        defaultHumanSoleus.curves.forceLengthDistalTitinInverseCurve,...                    
                        defaultHumanSoleus.curves.forceLengthIgPTitinCurve,...
                        defaultHumanSoleus.curves.forceLengthIgPTitinInverseCurve,...
                        defaultHumanSoleus.curves.forceLengthPevkTitinCurve,...
                        defaultHumanSoleus.curves.forceLengthPevkTitinInverseCurve,...
                        defaultHumanSoleus.curves.forceLengthIgDTitinCurve,...
                        defaultHumanSoleus.curves.forceLengthIgDTitinInverseCurve,...
                        defaultHumanSoleus.sarcomere.ZLineToT12NormLengthAtOptimalFiberLength,...
                        defaultHumanSoleus.sarcomere.IGPNormLengthAtOptimalFiberLength,...
                        defaultHumanSoleus.sarcomere.PEVKNormLengthAtOptimalFiberLength,...
                        defaultHumanSoleus.sarcomere.IGDFixedNormLengthAtOptimalFiberLength,...
                        defaultHumanSoleus.sarcomere.titinModelType);   
                      
      lengthZ2IgpHuman  = ...
          defaultHumanSoleus.sarcomere.ZLineToT12NormLengthAtOptimalFiberLength...
          .*ones(size(samples));
      lengthZ2PevkHuman     = lengthZ2IgpHuman + curveSampleIgpHuman.x;
      lengthZ2IgdHuman      = lengthZ2PevkHuman + curveSamplePevkHuman.x;  
      lengthZ2MyosinHuman   = lengthZ2IgdHuman + curveSampleIgdHuman.x;
      lengthCEHuman     = samples.*0.5;     

      fig_titinCurves=figure;
      
      loptHuman = defaultHumanSoleus.musculotendon.optimalFiberLength;
      colorIgp = [1,0,0];
      colorPEVK= [0.5,0,0];
      colorIgd = [1,0,0];

      plot( (lengthCEHuman.*2).*loptHuman,...
            (lengthZ2PevkHuman).*loptHuman,...
            '-','Color',colorIgp,...
            'LineWidth',1,...
            'DisplayName','Model: ZL to IgD/PEVK');
          
      hold on;
      
    
            
      plot( (lengthCEHuman.*2).*loptHuman,...
            (lengthZ2IgdHuman).*loptHuman,...
            '-','Color',colorPEVK,...
            'LineWidth',1,...
            'DisplayName','Model: ZL to PEVK/IgD');
          
            
      plot( (lengthCEHuman.*2).*loptHuman,...
            (lengthZ2MyosinHuman).*loptHuman,...
            '--','Color',colorIgd,...
            'LineWidth',1,...
            'DisplayName','Model: ZL to IgD/Myosin');

end





%%
% Convert all of the Bezier curves to quadratic splines Bezier curves and write the
% information to file
%%

felineSoleusNormMuscleQuadraticCurves=[];
humanSoleusNormMuscleQuadraticCurves=[];

curveNames = fieldnames(defaultFelineSoleus.curves);

indexCurve=1;
while indexCurve < length(curveNames)

    if(contains(curveNames{indexCurve},'use')==1 ...
        || isempty(defaultFelineSoleus.curves.(curveNames{indexCurve}))==1 ...
        || isstruct(defaultFelineSoleus.curves.(curveNames{indexCurve}))==0)
         curveNames(indexCurve)=[];
         indexCurve=indexCurve-1;
    end
        indexCurve=indexCurve+1;
end

fprintf('\n\nConverting 6th order curves to 4th order\n\n');
for indexCurve=1:1:length(curveNames)
        disp(curveNames{indexCurve});

        felineSoleusNormMuscleQuadraticCurves.(curveNames{indexCurve}) = ...
            convertToQuadraticBezierCurve(...
                defaultFelineSoleus.curves.(curveNames{indexCurve}),...
                numberOfQuadraticSubdivisions);

        humanSoleusNormMuscleQuadraticCurves.(curveNames{indexCurve}) = ...
            convertToQuadraticBezierCurve(...
                defaultHumanSoleus.curves.(curveNames{indexCurve}),...
                numberOfQuadraticSubdivisions);        
end



%%
% Convert all of the Bezier curves to quadratic Bezier curves and write the
% information to file
%%

humanSkeletalMuscleQuadraticCurves=[];
zeroHumanTitinQuadraticCurves=[];
oneHumanTitinQuadraticCurves=[];
curveNames = fieldnames(defaultHumanSoleus.curves);

zeroFelineTitinQuadraticCurves=[];
oneFelineTitinQuadraticCurves=[];


curveNamesTitin = [];
for i=1:1:length(curveNames)
    if(contains(curveNames{i},'Titin') && contains(curveNames{i},'use')==0)
        curveNamesTitin = [curveNamesTitin,{curveNames{i}}];
    end
end



indexCurve=1;
while indexCurve < length(curveNames)

    if(contains(curveNames{indexCurve},'use')==1 ...
        || isempty(defaultFelineSoleus.curves.(curveNames{indexCurve}))==1 ...
        || isstruct(defaultFelineSoleus.curves.(curveNames{indexCurve}))==0)
         curveNames(indexCurve)=[];
         indexCurve=indexCurve-1;
    end

    
        indexCurve=indexCurve+1;
end

fprintf('\n\nConverting 6th order curves to 4th order (human)\n\n');

for indexCurve=1:1:length(curveNames)
    disp(curveNames{indexCurve});
    humanSkeletalMuscleQuadraticCurves.(curveNames{indexCurve}) = ...
        convertToQuadraticBezierCurve(...
            defaultHumanSoleus.curves.(curveNames{indexCurve}),...
            numberOfQuadraticSubdivisions);
end


fprintf('\n\nConverting 6th order curves to 4th order (human-one-titin)\n\n');

for indexCurve=1:1:length(curveNamesTitin)
    disp(curveNamesTitin{indexCurve});
    oneHumanTitinQuadraticCurves.(curveNamesTitin{indexCurve}) = ...
        convertToQuadraticBezierCurve(...
            oneHumanSoleus.curves.(curveNamesTitin{indexCurve}),...
            numberOfQuadraticSubdivisions);
end

fprintf('\n\nConverting 6th order curves to 4th order (human-zero-titin)\n\n');

for indexCurve=1:1:length(curveNamesTitin)
    disp(curveNamesTitin{indexCurve});
    zeroHumanTitinQuadraticCurves.(curveNamesTitin{indexCurve}) = ...
        convertToQuadraticBezierCurve(...
            zeroHumanSoleus.curves.(curveNamesTitin{indexCurve}),...
            numberOfQuadraticSubdivisions);
end

fprintf('\n\nConverting 6th order curves to 4th order (feline-one-titin)\n\n');

for indexCurve=1:1:length(curveNamesTitin)
    disp(curveNamesTitin{indexCurve});
    oneFelineTitinQuadraticCurves.(curveNamesTitin{indexCurve}) = ...
        convertToQuadraticBezierCurve(...
            oneFelineSoleus.curves.(curveNamesTitin{indexCurve}),...
            numberOfQuadraticSubdivisions);
end

fprintf('\n\nConverting 6th order curves to 4th order (feline-zero-titin)\n\n');

for indexCurve=1:1:length(curveNamesTitin)
    disp(curveNamesTitin{indexCurve});
    zeroFelineTitinQuadraticCurves.(curveNamesTitin{indexCurve}) = ...
        convertToQuadraticBezierCurve(...
            zeroFelineSoleus.curves.(curveNamesTitin{indexCurve}),...
            numberOfQuadraticSubdivisions);
end

%Add the titin curves for interpolation

%%
% Plot everything
%%

colorOverride               = [1,0.5,0.5];
lineWidthOverride           = 1.5;
plotNumericalDerivatives    = 0;

if(flag_plotQuadraticFelineSoleus==1)
    fprintf('\n\nPlotting quadratic Bezier splines: feline soleus\n\n');

    structOfFigures = plotStructOfQuadraticBezierSplines( structOfFigures,...
                        felineSoleusNormMuscleQuadraticCurves,'Inverse',...
                        colorOverride,lineWidthOverride,...
                        'Feline soleus (2nd order)',2,...
                        plotNumericalDerivatives);
    
    colorOverride               = [0,0,0];
    lineWidthOverride           = 1;
    plotNumericalDerivatives    = 0;
end

if(flag_plotQuadraticHumanSoleus==1)
    fprintf('\n\nPlotting quadratic Bezier splines: human soleus\n\n');
    
    structOfFigures = plotStructOfQuadraticBezierSplines( structOfFigures,...
                        humanSkeletalMuscleQuadraticCurves,'Inverse',...
                        colorOverride,lineWidthOverride,...
                        'Human (2nd order)',3,...
                        plotNumericalDerivatives);
end


quadraticBezierFortranCurveFelineFolder = ['output/tables/FortranExport/QuadraticBezierFelineCurves/'];
quadraticBezierCsvCurveFelineFolder     = ['output/tables/curves/QuadraticBezierFelineCurves/'];


quadraticBezierFortranCurveHumanFolder = ['output/tables/FortranExport/QuadraticBezierHumanCurves/'];
quadraticBezierCsvCurveHumanFolder     = ['output/tables/curves/QuadraticBezierHumanCurves/'];


quadraticBezierFortranCurveFolderZero = ['output/tables/FortranExport/QuadraticBezierCurvesZero/'];
quadraticBezierCsvCurveFolderZero     = ['output/tables/curves/QuadraticBezierCurvesZero/'];

quadraticBezierFortranCurveFolderOne = ['output/tables/FortranExport/QuadraticBezierCurvesOne/'];
quadraticBezierCsvCurveFolderOne     = ['output/tables/curves/QuadraticBezierCurvesOne/'];



%%
% Save Matlab files
%%

save('output/structs/defaultFelineSoleusQuadraticCurves.mat', ...
    'felineSoleusNormMuscleQuadraticCurves');

save('output/structs/defaultHumanSoleusQuadraticCurves.mat', ...
    'humanSoleusNormMuscleQuadraticCurves');

%We also output a version of the human curves assuming that the proximal
%end of the PEVK segment attaches to actin, and another where the distal
%end of the PEVK segment attaches to actin. The coefficients for these two
%curves are interpolated to form, on the fly, the curve that defines the
%force-length characteristics of titin given an attachement point loaction
%anywhere within the PEVK segment.

%Contains titin curves assuming that the proximal end of the PEVK 
%segment attaches to actin
save('output/structs/defaultHumanSoleusQuadraticCurves_TitinActinPEVKProxAttach.mat', ...
    'zeroHumanTitinQuadraticCurves');

%Contains titin curves assuming that the distal end of the PEVK 
%segment attaches to actin
save('output/structs/defaultHumanSoleusQuadraticCurves_TitinActinPEVKDistalAttach.mat', ...
    'oneHumanTitinQuadraticCurves');

%%
% Generate the Fortran code to create each curve ...
%%

fprintf('\n\nwriting Fortran code: feline soleus\n\n');

status = writeBezierCurveStructuresToFortran(...
    felineSoleusNormMuscleQuadraticCurves,...
    'felineSoleus',...
    quadraticBezierFortranCurveFelineFolder);

fprintf('\n\nwriting Fortran code: human soleus\n\n');

status = writeBezierCurveStructuresToFortran(...
    humanSoleusNormMuscleQuadraticCurves,...
    'humanSkeletalMuscle',...
    quadraticBezierFortranCurveHumanFolder);


fprintf('\n\nwriting Fortran code: human soleus titin (one)\n\n');

status = writeBezierCurveStructuresToFortran(...
    oneHumanTitinQuadraticCurves,...
    'humanSkeletalMuscleTitinOne',...
    quadraticBezierFortranCurveFolderOne);

fprintf('\n\nwriting Fortran code: human soleus titin (zero)\n\n');

status = writeBezierCurveStructuresToFortran(...
    zeroHumanTitinQuadraticCurves,...
    'humanSkeletalMuscleTitinZero',...
    quadraticBezierFortranCurveFolderZero);

fprintf('\n\nwriting Fortran code: feline soleus titin (one)\n\n');

status = writeBezierCurveStructuresToFortran(...
    oneFelineTitinQuadraticCurves,...
    'felineSkeletalMuscleTitinOne',...
    quadraticBezierFortranCurveFolderOne);

fprintf('\n\nwriting Fortran code: feline soleus titin (zero)\n\n');

status = writeBezierCurveStructuresToFortran(...
    zeroFelineTitinQuadraticCurves,...
    'felineSkeletalMuscleTitinZero',...
    quadraticBezierFortranCurveFolderZero);

%%
% CSV file of the control points for each curve ...
%%

fprintf('\n\nwriting Bezier curves to CSV: feline soleus\n\n');

status = writeBezierCurveStructuresToCSV(...
    felineSoleusNormMuscleQuadraticCurves,...
    'felineSoleus',...
    quadraticBezierCsvCurveFelineFolder);


fprintf('\n\nwriting Bezier curves to CSV: human soleus\n\n');

status = writeBezierCurveStructuresToCSV(...
    humanSkeletalMuscleQuadraticCurves,...
    'humanSkeletalMuscle',...
    quadraticBezierCsvCurveHumanFolder);  

muscleStruct = readBezierCurveStructuresFromCSV(quadraticBezierCsvCurveHumanFolder);

fprintf('\n\nwriting Bezier curves to CSV: human titin soleus (zero)\n\n');

status = writeBezierCurveStructuresToCSV(...
    zeroHumanTitinQuadraticCurves,...
    'humanSkeletalMuscleTitinZero',...
    quadraticBezierCsvCurveFolderZero);  

fprintf('\n\nwriting Bezier curves to CSV: human titin soleus (one)\n\n');

status = writeBezierCurveStructuresToCSV(...
    oneHumanTitinQuadraticCurves,...
    'humanSkeletalMuscleTitinOne',...
    quadraticBezierCsvCurveFolderOne);  

fprintf('\n\nwriting Bezier curves to CSV: feline titin soleus (zero)\n\n');

status = writeBezierCurveStructuresToCSV(...
    zeroFelineTitinQuadraticCurves,...
    'felineSkeletalMuscleTitinZero',...
    quadraticBezierCsvCurveFolderZero);  

fprintf('\n\nwriting Bezier curves to CSV: feline titin soleus (one)\n\n');

status = writeBezierCurveStructuresToCSV(...
    oneFelineTitinQuadraticCurves,...
    'felineSkeletalMuscleTitinOne',...
    quadraticBezierCsvCurveFolderOne);  

%%
% Generate a numerical sample of each curve
%%

quadraticBezierCurveFelineSamplesFolder = ['output/sample/QuadraticBezierFelineCurves/'];
quadraticBezierCurveHumanSamplesFolder  = ['output/sample/QuadraticBezierHumanCurves/'];

nPts=100;


status = writeQuadraticBezierCurveSamplesToCSV(...
            felineSoleusNormMuscleQuadraticCurves,...
            nPts,...
            'felineSoleus',...
            quadraticBezierCurveFelineSamplesFolder);

status = writeQuadraticBezierCurveSamplesToCSV(...
            zeroFelineTitinQuadraticCurves,...
            nPts,...
            'felineSoleusZero',...
            quadraticBezierCurveFelineSamplesFolder);

status = writeQuadraticBezierCurveSamplesToCSV(...
            oneFelineTitinQuadraticCurves,...
            nPts,...
            'felineSoleusOne',...
            quadraticBezierCurveFelineSamplesFolder);



status = writeQuadraticBezierCurveSamplesToCSV(...
            humanSoleusNormMuscleQuadraticCurves,...
            nPts,...
            'humanSoleus',...
            quadraticBezierCurveHumanSamplesFolder);

status = writeQuadraticBezierCurveSamplesToCSV(...
            zeroHumanTitinQuadraticCurves,...
            nPts,...
            'humanSoleusZero',...
            quadraticBezierCurveHumanSamplesFolder);

status = writeQuadraticBezierCurveSamplesToCSV(...
            oneHumanTitinQuadraticCurves,...
            nPts,...
            'humanSoleusOne',...
            quadraticBezierCurveHumanSamplesFolder);



%%
% Remove the directories ...
%%
rmpath(parametersDirectoryTreeMTParams);
rmpath(parametersDirectoryTreeExperiments);
rmpath(parametersDirectoryTreeModels);
rmpath(parametersDirectoryTreeCurves);
