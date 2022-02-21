

clc;
close all;
clear all;
fitCrossBridgeStiffnessDampingToKirch199490Hz=0;
flag_useFixedLambdaECM    = 0;

flag_plotEveryCurve       = 0;
flag_plotForceLengthDetail= 0;


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

if(flag_plotEveryCurve==1)
    figH = plotStructOfBezierSplines( felineSoleusNormMuscleCurves,...
                                      'Inverse');                          
    
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
% Convert all of the Bezier curves to hyperbolic splines and write the
% information to file
%%

n =size(felineSoleusNormMuscleCurves.activeForceLengthCurve.xpts,2);

activeForceLengthCurveH = struct(...
    'xpts',zeros(2,n),...
    'ypts',zeros(4,n),...
    'xEnd',felineSoleusNormMuscleCurves.activeForceLengthCurve.xEnd,...
    'yEnd',felineSoleusNormMuscleCurves.activeForceLengthCurve.yEnd,...
    'dydxEnd',felineSoleusNormMuscleCurves.activeForceLengthCurve.dydxEnd,...
    'd2ydx2End',zeros(1,2).*NaN);

xRange = felineSoleusNormMuscleCurves.activeForceLengthCurve.xEnd;

for i=1:1:n

    x0 = felineSoleusNormMuscleCurves.activeForceLengthCurve.xpts(1,i);
    x1 = felineSoleusNormMuscleCurves.activeForceLengthCurve.xpts(end,i);
    
    y0 = felineSoleusNormMuscleCurves.activeForceLengthCurve.ypts(1,i);
    y1 = felineSoleusNormMuscleCurves.activeForceLengthCurve.ypts(end,i);

    dydx0 = calcBezierYFcnXDerivative(x0,...
        felineSoleusNormMuscleCurves.activeForceLengthCurve,1);

    dydx1 = calcBezierYFcnXDerivative(x1,...
        felineSoleusNormMuscleCurves.activeForceLengthCurve,1);

    dudx = 1/(x1-x0);
    dxdu = 1/dudx;

    dydu0 = dydx0*dxdu;
    dydu1 = dydx1*dxdu;



    d = 1;

    sqrt_dydu0dydu1=sqrt(dydu0*dudu1);

    c0 = -(sqrt(dydu0*dydu1)+dydu1)/dydu1;
    c1 =  (sqrt(dydu0*dydu1)-dydu1)/dydu1;
    c  = c0;
    b = (sqrt_dydu0dydu1*y1 + dydu*dydu0)/sqrt_dydu0dydu1;

    a0 = dydu0+b*c0;
    a1 = dydu0+b*c1;
    a  = a0;

    err0 = (b0+a0)/(c+1) - y1;
    err1 = (b1+a1)/(c+1) - y1;
    if(abs(err1) < abs(err0))        
        a = a1;
        c = c1;
    end


    %Check
    u =0;
    y   = (a*u+b)/(c*u+1);
    t0  = (c*u+1);
    dudu= (a/t0)-(c*(a*u+b)/(t0*t0);
    
    assert(abs(y1c-y1)       < eps*10);
    assert(abs(dydu1c-dydu1) < eps*10);

    %Update the structure
    activeForceLengthCurveH.xpts(:,n)=[x0;x1];
    activeForceLengthCurveH.ypts(:,n)=[a;b;c;d];
    if(i==1)
        activeForceLengthCurveH.xEnd(1,1)=x0;        
        activeForceLengthCurveH.yEnd(1,1)=y0;
        activeForceLengthCurveH.dydxEnd(1,1)=dydu0*dudx;
        activeForceLengthCurveH.d2ydx2End(1,1)=d2ydu20c*(dudx*dudx);
    end
    if(i == n)
        activeForceLengthCurveH.xEnd(1,2)=x1;        
        activeForceLengthCurveH.yEnd(1,2)=y1;
        activeForceLengthCurveH.dydxEnd(1,2)=dydu1*dudx;
        activeForceLengthCurveH.d2ydx2End(1,2)=d2ydu21c*(dudx*dudx);
    end
end

npts=100;
xSample= [0:(1/(npts-1)):1].*(xRange(1,2)-xRange(1,1)) + xRange(1,1);
ySample = zeros(npts,3);


%felineSoleusNormMuscleHyperbolicCurves =felineSoleusNormMuscleCurves;

% for i=1:1:length(curveNames)
%     if(isempty(felineSoleusNormMuscleCurves.(curveNames{i}))==0)
%         numberOfSegments = size(felineSoleusNormMuscleCurves.(curveNames{i}).xpts,2);
%         felineSoleusNormMuscleHyperbolicCurves.(curveNames{i}) = ...
%             struct('a', size())
% 
%      = createHyperbolicSplineFromBezierSpline(...
%          felineSoleusNormMuscleCurves.(curveNames{i}));
%     end
% end

   
%%
% Remove the directories ...
%%
rmpath(parametersDirectoryTreeMTParams);
rmpath(parametersDirectoryTreeExperiments);
rmpath(parametersDirectoryTreeModels);
rmpath(parametersDirectoryTreeCurves);
