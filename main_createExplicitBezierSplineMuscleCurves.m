

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



scaleOptimalFiberLength      = 1.0; 

scaleMaximumIsometricTension = 1;

if(exist('fitCrossBridgeStiffnessDampingToKirch199490Hz','var')==0)
  fitCrossBridgeStiffnessDampingToKirch199490Hz = 1;
end

%%
% Fitted properties for a feline soleus
%%

[felineSoleusMusculotendonProperties, ...
 felineSoleusSarcomereProperties,...
 felineSoleusActiveForceLengthData,...
 felineSoleusPassiveForceLengthData] = createFelineSoleus(...                                          
                                          scaleOptimalFiberLength,...
                                          scaleMaximumIsometricTension,...
                                          fitCrossBridgeStiffnessDampingToKirch199490Hz,...
                                          flag_useOctave);

createMusculoTendonFcn = ...
  @(argScaleFiberLength,argScaleFiso)createFelineSoleus(...
                                        argScaleFiberLength,...
                                        argScaleFiso,...
                                        flag_useOctave); 
      
disp('Set to match the force-velocity curve of umat41 (EHTMM)');
felineSoleusMusculotendonProperties.forceVelocityMultiplierAtHalfMaximumFiberVelocity=0.25;

[felineSoleusNormMuscleCurves,...
 felineSoleusMusculotendonPropertiesUpd,...
 felineSoleusSarcomerePropertiesUpd,...
 activeForceLengthCurveAnnotationPoints,...
 felineSoleusActiveForceLengthDataUpd,...
 felineSoleusPassiveForceLengthDataUpd,...
 forceLengthCurveSettings]= ...
    createFittedMuscleCurves( ...
      felineSoleusMusculotendonProperties,...
      felineSoleusSarcomereProperties,...
      felineSoleusActiveForceLengthData,...
      felineSoleusPassiveForceLengthData,...
      shiftLengthActiveForceLengthCurveDescendingCurve,...
      flag_useFixedLambdaECM,...
      flag_enableNumericallyNonZeroGradients,...
      smallNumericallyNonZeroNumber,...
      flag_solveForOptimalFiberLengthOfBestFit,...
      createMusculoTendonFcn,...
      flag_useOctave);

lambdaECM = felineSoleusSarcomerePropertiesUpd.extraCellularMatrixPassiveForceFraction;    

curveNames = fields(felineSoleusNormMuscleCurves);
for indexCurve = 1:1:length(curveNames)
    if(isstruct(felineSoleusNormMuscleCurves.(curveNames{indexCurve}))==1)
        if(isempty(felineSoleusNormMuscleCurves.(curveNames{indexCurve}))==0)
            if(contains( curveNames{indexCurve},'ECM'))
                felineSoleusNormMuscleCurves.(curveNames{indexCurve}) = ...
                    scaleCurveStruct(1,(1/lambdaECM),...
                        felineSoleusNormMuscleCurves.(curveNames{indexCurve}));
            end
            if(contains( curveNames{indexCurve},'Titin'))
                felineSoleusNormMuscleCurves.(curveNames{indexCurve}) = ...
                    scaleCurveStruct(1,(1/(1-lambdaECM)),...
                       felineSoleusNormMuscleCurves.(curveNames{indexCurve}));
            end
            
        end
    end
end


humanSkeletalMuscleCurves =struct('activeForceLengthCurve',[]);
%Get the default sarcomere properties for a feline soles          
scaleOptimalFiberLength                         = 1;
fitCrossBridgeStiffnessDampingToKirch199490Hz   = 1;
flag_Cat1_Human2                                = 2;                           

[humanSoleusSarcomereProperties] =...
  getMammalianSkeletalMuscleNormalizedSarcomereProperties(...
    scaleOptimalFiberLength,...
    flag_Cat1_Human2,...
    fitCrossBridgeStiffnessDampingToKirch199490Hz);

curvinessActiveForceLength              = 1;
flag_compensateForCrossbridgeStiffness  = 1;
computeIntegral=0;
[humanSkeletalMuscleCurves.activeForceLengthCurve, ...
  activeForceLengthCurveAnnotationPoints] ...
    = createFiberActiveForceLengthCurve(...
          humanSoleusSarcomereProperties.normMyosinHalfLength*2,...
          humanSoleusSarcomereProperties.normMyosinBareHalfLength*2,...
          humanSoleusSarcomereProperties.normActinLength,...
          humanSoleusSarcomereProperties.normZLineLength,...
          humanSoleusSarcomereProperties.normSarcomereLengthZeroForce,...
          humanSoleusSarcomereProperties.normCrossBridgeStiffness,...
          curvinessActiveForceLength, ... 
          shiftLengthActiveForceLengthCurveDescendingCurve,...
          flag_compensateForCrossbridgeStiffness,...
          flag_enableNumericallyNonZeroGradients,...
          smallNumericallyNonZeroNumber,...
          computeIntegral,...
          'humanSkeletalMuscle',...
          flag_useOctave); 

%%
%Plot feline curves + experimental data from Herzog & Leonard 2002
%%

structOfFigures = [];

if(flag_plotEveryCurve==1)
    figH = plotStructOfBezierSplines( structOfFigures,...
                                      felineSoleusNormMuscleCurves,...
                                      {'Inverse','use'},[],[],[],[],0);                          
    
    %%
    % Note the average offset between the active-force-length curve and
    % the transformed data
    %%
    
    xExp = felineSoleusActiveForceLengthDataUpd(2:end,1);
    yExp = felineSoleusActiveForceLengthDataUpd(2:end,2);
    xCurve = zeros(size(xExp));
    
    for i=1:1:length(xExp)
    xCurve(i,1) = calcBezierFcnXGivenY(yExp(i,1), ...
      felineSoleusNormMuscleCurves.activeForceLengthCurve,... 
      xExp(i,1));
    end                                    
    dx = mean(xCurve-xExp);
    %felineSoleusActiveForceLengthDataUpd(:,1)=...
    %  felineSoleusActiveForceLengthDataUpd(:,1)+dx;
    
    disp('Normalized length offset');
    fprintf('%1.6f lce/lopt \n',felineSoleusActiveForceLengthDataUpd(1,1));
    disp('Average error on the descending limb');
    fprintf('%1.6f lce/lopt \n',dx);
    fprintf('%1.6f mm \n',dx*(felineSoleusMusculotendonPropertiesUpd.optimalFiberLength*1000));
    
    lceNStart = felineSoleusActiveForceLengthDataUpd(1,1);
    save('output/structs/normalizedFiberLengthStartHerzogLeonard2002.mat',...
         'lceNStart');
    
    % dl = felineSoleusActiveForceLengthDataUpd(2:end,1)-1;
    % A  = [dl ones(size(dl))];
    % b  = felineSoleusActiveForceLengthDataUpd(2:end,2);
    % 
    % x     = (A'*A)\(A'*b);
    % y0    = x(2,1);
    % dydx0 = x(1,1);
    % 
    % felineSoleusActiveForceLengthLineBestFit = zeros(length(dl),1);
    % 
    % felineSoleusActiveForceLengthLineBestFit(:,1) = ...
    %   felineSoleusActiveForceLengthDataUpd(2:end,1);
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
                        felineSoleusNormMuscleCurves.('tendonForceLengthCurve'), 200,[]);
    
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
      plot(  felineSoleusActiveForceLengthDataUpd(:,1),...
           felineSoleusActiveForceLengthDataUpd(:,2),'xb');
      hold on;          
    
    figure(figH.fiberForceLengthCurve);
      subplot(2,2,1);
      plot(   felineSoleusPassiveForceLengthDataUpd(:,1),...
              felineSoleusPassiveForceLengthDataUpd(:,2),'xb');
      hold on;          
  
end

%%
% Plot the passive force length curve of the 2 segment titin model and
% compare it to the passive force length curve
%%
if(flag_plotForceLengthDetail==1)
    lceN0 = calcBezierFcnXGivenY(0, ...
              felineSoleusNormMuscleCurves.fiberForceLengthCurve);
    lceN1 = calcBezierFcnXGivenY(1, ...
              felineSoleusNormMuscleCurves.fiberForceLengthCurve);
            
    npts = 100;
    lceNSeries = [lceN0:((lceN1-lceN0)/(npts-1)):lceN1]';
    
    fNSeries = zeros(npts,4); % fpe, fecm, fIgp, fPevkIgd,
    lNSeries = zeros(npts,4); % lpe, lecm, lIgp, lPevkIgd,
    
    normLengthIgdFixed = ...
      felineSoleusSarcomerePropertiesUpd.IGDFixedNormLengthAtOptimalFiberLength;
    
    normLengthT12ToZ = ...
      felineSoleusSarcomerePropertiesUpd.ZLineToT12NormLengthAtOptimalFiberLength;
    
    
    for i=1:1:npts
      lceN  = lceNSeries(i,1);
      fpeN  = calcBezierYFcnXDerivative(lceN,...
                felineSoleusNormMuscleCurves.fiberForceLengthCurve,0);
      fecmN = calcBezierYFcnXDerivative(lceN*0.5,...
                felineSoleusNormMuscleCurves.forceLengthECMHalfCurve,0);
              
      lIgpPevkN = lceN*0.5 - normLengthIgdFixed - normLengthT12ToZ; 
              
      [lIgpN, lPevkIgdN, fTiN] = calcSeriesSpringStretch(lIgpPevkN,...
                felineSoleusNormMuscleCurves.forceLengthIgpCurve,...
                felineSoleusNormMuscleCurves.forceLengthIgpInverseCurve, ...
                felineSoleusNormMuscleCurves.forceLengthPevkIgdCurve,...
                felineSoleusNormMuscleCurves.forceLengthPevkIgdInverseCurve);          
    
      fNSeries(i,:) = [fpeN,fecmN,fTiN, fTiN];
      lNSeries(i,:) = [lceN, lceN,lIgpN*2, lPevkIgdN*2];
              
    end
    
    
    fig_forceLength = figure;
      plot(lNSeries(:,1),fNSeries(:,1),'k');
      hold on;
      plot(lNSeries(:,1), fNSeries(:,2)+fNSeries(:,3),'b');
      hold on;
      plot(lNSeries(:,1), fNSeries(:,2),'r');
      hold on;
      plot(felineSoleusPassiveForceLengthDataUpd(:,1),...
           felineSoleusPassiveForceLengthDataUpd(:,2),'xb');
      legend('fpe','fecm+fti','fecm','data');
      xlabel('Norm. Length')
      ylabel('Norm. Force');
      

end

defaultFelineSoleus = struct('musculotendon',...
                            felineSoleusMusculotendonPropertiesUpd,...
                            'sarcomere',...
                            felineSoleusSarcomerePropertiesUpd,...
                            'falData',...
                            felineSoleusActiveForceLengthDataUpd,...
                            'fpeData',...
                            felineSoleusPassiveForceLengthDataUpd,...
                            'curves',...
                            felineSoleusNormMuscleCurves);
                      
save('output/structs/defaultFelineSoleus.mat',...
     'defaultFelineSoleus');                      



%%
% Convert all of the Bezier curves to quadratic splines Bezier curves and write the
% information to file
%%

felineSoleusNormMuscleQuadraticCurves=[];
curveNames = fieldnames(felineSoleusNormMuscleCurves);

indexCurve=1;
while indexCurve < length(curveNames)

    if(contains(curveNames{indexCurve},'use')==1 ...
        || isempty(felineSoleusNormMuscleCurves.(curveNames{indexCurve}))==1 ...
        || isstruct(felineSoleusNormMuscleCurves.(curveNames{indexCurve}))==0)
         curveNames(indexCurve)=[];
         indexCurve=indexCurve-1;
    end
        indexCurve=indexCurve+1;
end


for indexCurve=1:1:length(curveNames)
        disp(curveNames{indexCurve});

        felineSoleusNormMuscleQuadraticCurves.(curveNames{indexCurve}) = ...
            convertToQuadraticBezierCurve(...
                felineSoleusNormMuscleCurves.(curveNames{indexCurve}),...
                numberOfQuadraticSubdivisions);
    
end


%%
% Convert all of the Bezier curves to cubic splines Bezier curves and write the
% information to file
%%
% felineSoleusNormMuscleCubicCurves=[];
% curveNames = fieldnames(felineSoleusNormMuscleCurves);
% 
% indexCurve=1;
% while indexCurve < length(curveNames)
% 
%     if(contains(curveNames{indexCurve},'use')==1 ...
%         || isempty(felineSoleusNormMuscleCurves.(curveNames{indexCurve}))==1 ...
%         || isstruct(felineSoleusNormMuscleCurves.(curveNames{indexCurve}))==0)
%          curveNames(indexCurve)=[];
%          indexCurve=indexCurve-1;
%     end
%         indexCurve=indexCurve+1;
% end

% for indexCurve=1:1:length(curveNames)
%     disp(curveNames{indexCurve});
%     felineSoleusNormMuscleCubicCurves.(curveNames{indexCurve}) = ...
%         convertToCubicBezierCurve(...
%             felineSoleusNormMuscleCurves.(curveNames{indexCurve}),...
%             numberOfCubicSubdivisions,cubicZeroSecondDerivative,1);
% end


%%
% Convert all of the Bezier curves to quadratic Bezier curves and write the
% information to file
%%

humanSkeletalMuscleQuadraticCurves=[];
curveNames = fieldnames(humanSkeletalMuscleCurves);

indexCurve=1;
while indexCurve < length(curveNames)

    if(contains(curveNames{indexCurve},'use')==1 ...
        || isempty(felineSoleusNormMuscleCurves.(curveNames{indexCurve}))==1 ...
        || isstruct(felineSoleusNormMuscleCurves.(curveNames{indexCurve}))==0)
         curveNames(indexCurve)=[];
         indexCurve=indexCurve-1;
    end
        indexCurve=indexCurve+1;
end

for indexCurve=1:1:length(curveNames)
    disp(curveNames{indexCurve});
    humanSkeletalMuscleQuadraticCurves.(curveNames{indexCurve}) = ...
        convertToQuadraticBezierCurve(...
            humanSkeletalMuscleCurves.(curveNames{indexCurve}),...
            numberOfQuadraticSubdivisions);
end

%%
% Plot everything
%%

%structOfFigures = [];
% colorOverride               = [1,1,1].*0.75;
% lineWidthOverride           = 2;
% plotNumericalDerivatives    = 0;
% structOfFigures = plotStructOfBezierSplines( structOfFigures, ...
%                     felineSoleusNormMuscleCurves,'Inverse',...
%                     colorOverride,lineWidthOverride,...
%                     'Feline soleus (5th order)',0,...
%                     plotNumericalDerivatives);
% 
% colorOverride               = [0.5,0.5,1];
% lineWidthOverride           = 1.5;
% plotNumericalDerivatives    = 0;
% 
% structOfFigures = plotStructOfBezierSplines( structOfFigures,...
%                     felineSoleusNormMuscleCubicCurves,'Inverse',...
%                     colorOverride,lineWidthOverride,...
%                     'Feline soleus (3rd order)',1,...
%                     plotNumericalDerivatives);
% 
colorOverride               = [1,0.5,0.5];
lineWidthOverride           = 1.5;
plotNumericalDerivatives    = 0;

structOfFigures = plotStructOfQuadraticBezierSplines( structOfFigures,...
                    felineSoleusNormMuscleQuadraticCurves,'Inverse',...
                    colorOverride,lineWidthOverride,...
                    'Feline soleus (2nd order)',2,...
                    plotNumericalDerivatives);

colorOverride               = [0,0,0];
lineWidthOverride           = 1;
plotNumericalDerivatives    = 0;
structOfFigures = plotStructOfQuadraticBezierSplines( structOfFigures,...
                    humanSkeletalMuscleQuadraticCurves,'Inverse',...
                    colorOverride,lineWidthOverride,...
                    'Human (2nd order)',3,...
                    plotNumericalDerivatives);


%Write the fortran matrices
quadraticBezierFortranCurveFolder = ['output/tables/QuadraticBezierFortranCurves/'];
quadraticBezierCsvCurveFolder     = ['output/tables/QuadraticBezierCSVCurves/'];

status = writeMuscleStructuresToFortran(...
    felineSoleusNormMuscleQuadraticCurves,...
    'felineSoleus',...
    quadraticBezierFortranCurveFolder);

status = writeMuscleStructuresToFortran(...
    humanSkeletalMuscleCurves,...
    'humanSkeletalMuscle',...
    quadraticBezierFortranCurveFolder);

status = writeMuscleStructuresToCSV(...
    felineSoleusNormMuscleQuadraticCurves,...
    'felineSoleus',...
    quadraticBezierCsvCurveFolder);

status = writeMuscleStructuresToCSV(...
    humanSkeletalMuscleQuadraticCurves,...
    'humanSkeletalMuscle',...
    quadraticBezierCsvCurveFolder);   

muscleStruct = readMuscleStructuresFromCSV(quadraticBezierCsvCurveFolder);


%%
% Remove the directories ...
%%
rmpath(parametersDirectoryTreeMTParams);
rmpath(parametersDirectoryTreeExperiments);
rmpath(parametersDirectoryTreeModels);
rmpath(parametersDirectoryTreeCurves);
