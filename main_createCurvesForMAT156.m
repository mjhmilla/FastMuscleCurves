

clc;
close all;
clear all;

opengl('save','software');

%Warning: This option is not working properly. The code will run,
%         curves will be produced, but when used in simulation the 
%         resulting curves do not behave as expected.
flag_adjustCurvesImplicitlyIncludeTendon = 0;

flag_solveForHerzogLeonard1997Params=1;
flag_addTendonStretchtToFpe=0;

%%
% Architectural Parameters from the LS-DYNA files
%%

catSoleusHL2002.lceOpt  =0.0428571;
catSoleusHL2002.fceOpt  =21.612138;
catSoleusHL2002.lceOptAT=0.0425376;
catSoleusHL2002.fceOptAT=21.451044;
catSoleusHL2002.lmtOptAT=0.0729888;
catSoleusHL2002.penOpt  =0.1221730;
catSoleusHL2002.penOptD =7.0 ;
catSoleusHL2002.ltSlk   =0.0304511;
catSoleusHL2002.et      =0.0458333;
catSoleusHL2002.vceMax  =4.50;

catSoleusHL1997.lceOpt  =0.038;
catSoleusHL1997.fceOpt  =41.661837;
catSoleusHL1997.lceOptAT=0.0377168;
catSoleusHL1997.fceOptAT=41.351296;
catSoleusHL1997.lmtOptAT=0.0681679;
catSoleusHL1997.penOpt  =0.1221730;
catSoleusHL1997.penOptD =7.0 ;
catSoleusHL1997.ltSlk   =0.0304511;
catSoleusHL1997.et      =0.0458333;
catSoleusHL1997.vceMax  =4.50;

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
parametersDirectoryTreePostProcessing= genpath('postprocessing');

addpath(parametersDirectoryTreeMTParams);
addpath(parametersDirectoryTreeExperiments);
addpath(parametersDirectoryTreeModels);
addpath(parametersDirectoryTreeCurves);
addpath(parametersDirectoryTreePostProcessing);

%%
% Load the curves that we've already generated for the cat soleus
%%

load('output/structs/defaultFelineSoleusQuadraticCurves.mat');
load('output/structs/defaultFelineSoleus.mat');
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

% Sample the fpe curve created by the ECM and series connection of titin
samples = fpeValues.x;

[ curveSampleECMHalfFeline,...
  curveSampleTitinFeline,...
  curveSampleTitinActiveFeline,...
  curveSampleIgpFeline,...
  curveSamplePevkFeline,...
  curveSampleIgdFeline,...
  curveSampleProximalTitinFeline,...
  curveSampleDistalTitinFeline] = ...
  sampleQuadraticTitinCurves(samples.*0.5,...
                    felineSoleusNormMuscleQuadraticCurves.forceLengthECMHalfCurve,...
                    felineSoleusNormMuscleQuadraticCurves.forceLengthProximalTitinCurve,...
                    felineSoleusNormMuscleQuadraticCurves.forceLengthProximalTitinInverseCurve,...
                    felineSoleusNormMuscleQuadraticCurves.forceLengthDistalTitinCurve,...
                    felineSoleusNormMuscleQuadraticCurves.forceLengthDistalTitinInverseCurve,...                    
                    felineSoleusNormMuscleQuadraticCurves.forceLengthIgPTitinCurve,...
                    felineSoleusNormMuscleQuadraticCurves.forceLengthIgPTitinInverseCurve,...
                    felineSoleusNormMuscleQuadraticCurves.forceLengthPevkTitinCurve,...
                    felineSoleusNormMuscleQuadraticCurves.forceLengthPevkTitinInverseCurve,...
                    felineSoleusNormMuscleQuadraticCurves.forceLengthIgDTitinCurve,...
                    felineSoleusNormMuscleQuadraticCurves.forceLengthIgDTitinInverseCurve,...
                    defaultFelineSoleus.sarcomere.ZLineToT12NormLengthAtOptimalFiberLength,...
                    defaultFelineSoleus.sarcomere.IGPNormLengthAtOptimalFiberLength,...
                    defaultFelineSoleus.sarcomere.PEVKNormLengthAtOptimalFiberLength,...
                    defaultFelineSoleus.sarcomere.IGDFixedNormLengthAtOptimalFiberLength,...
                    defaultFelineSoleus.sarcomere.titinModelType); 

fpe31Values.x = samples;
fpe31Values.y = curveSampleECMHalfFeline.y.*(defaultFelineSoleus.sarcomere.extraCellularMatrixPassiveForceFraction) ...
               +curveSampleTitinFeline.y.*(1-defaultFelineSoleus.sarcomere.extraCellularMatrixPassiveForceFraction);

%%
% Generate the values along the tendon
%%
fiso     = defaultFelineSoleus.musculotendon.fiso;
lceOpt   = defaultFelineSoleus.musculotendon.optimalFiberLength;
alphaOpt = defaultFelineSoleus.musculotendon.pennationAngle;
vMax     = lceOpt*defaultFelineSoleus.musculotendon.maximumNormalizedFiberVelocity;

lceOptAT = lceOpt*cos(alphaOpt);
ltSlk = defaultFelineSoleus.musculotendon.tendonSlackLength;

%lce*sin(alpha) = h
%dlce*sin(alpha)+lce*cos(alpha)*alphaDot=0
%alphaDot= -(dlce/lce)*tan(alpha)
alphaDotOpt = -(vMax/lceOpt)*tan(alphaOpt);
vMaxAT      = vMax*cos(alphaOpt)-lceOpt*sin(alphaOpt)*alphaDotOpt;
fisoAT      = fiso*cos(alphaOpt);


falValues.xAT   = zeros(size(falValues.x));
falValues.yAT   = zeros(size(falValues.y));

fpeValues.xAT   = zeros(size(fpeValues.x));
fpeValues.yAT   = zeros(size(fpeValues.y));

fpeAdjustedValues = fpeValues;
fpeAdjustedValues.xAT   = zeros(size(fpeValues.x));
fpeAdjustedValues.yAT   = zeros(size(fpeValues.y));

fvValues.xAT    = zeros(size(fvValues.x));
fvValues.yAT    = zeros(size(fvValues.y));

fpe31Values.xAT =zeros(size(fpe31Values.x));
fpe31Values.yAT =zeros(size(fpe31Values.y));


%Evaluate the active-force-length curve along the tendon
for i=1:1:length(falValues.x)

    if(falValues.x(i,1) >= 1)
        here=1;
    end

    fiberKinematics = calcFixedWidthPennatedFiberKinematicsAlongTendon(...
                                                falValues.x(i,1).*lceOpt,...
                                                0,...
                                                lceOpt,...
                                                alphaOpt);
    alpha = fiberKinematics.pennationAngle;
    lceAT = fiberKinematics.fiberLengthAlongTendon;

    falValues.xAT(i,1) = lceAT/lceOptAT;
    falValues.yAT(i,1) = fiso*falValues.y(i,1)*cos(alpha)/fisoAT;

    if(fiberKinematics.isClamped)
        fprintf('%i. fal isClamped \n',i);
    end
end

%Evaluate the passive-force-length curve along the tendon
for i=1:1:length(fpeValues.x)
    fiberKinematics = calcFixedWidthPennatedFiberKinematicsAlongTendon(...
                                                fpeValues.x(i,1).*lceOpt,...
                                                0,...
                                                lceOpt,...
                                                alphaOpt);

    alpha = fiberKinematics.pennationAngle;
    lceAT = fiberKinematics.fiberLengthAlongTendon;

    fpeValues.xAT(i,1) = lceAT/lceOptAT;
    fpeValues.yAT(i,1) = fiso*fpeValues.y(i,1)*cos(alpha)/fisoAT;

    if(fiberKinematics.isClamped)
        fprintf('%i. fpe isClamped \n',i);
    end    
end

%Evaluate the passive-force-length curve along the tendon
for i=1:1:length(fpeAdjustedValues.x)
    fiberKinematics = calcFixedWidthPennatedFiberKinematicsAlongTendon(...
                                                fpeAdjustedValues.x(i,1).*lceOpt,...
                                                0,...
                                                lceOpt,...
                                                alphaOpt);

    alpha = fiberKinematics.pennationAngle;
    lceAT = fiberKinematics.fiberLengthAlongTendon;
    
    fpeNAT = fiso*fpeValues.y(i,1)*cos(alpha)/fisoAT;

    ltN = calcQuadraticBezierYFcnXDerivative(fpeNAT,...
       felineSoleusNormMuscleQuadraticCurves.tendonForceLengthInverseCurve,0);
    ltDelta = ltN*ltSlk - ltSlk;

    if(flag_addTendonStretchtToFpe)
        fpeAdjustedValues.xAT(i,1) = (lceAT+ltDelta)/lceOptAT;
    else
        fpeAdjustedValues.xAT(i,1) = (lceAT)/lceOptAT;
    end
    fpeAdjustedValues.yAT(i,1) = fpeNAT;

    if(fiberKinematics.isClamped)
        fprintf('%i. fpe isClamped \n',i);
    end    
end

kpeWithTendonDeltaValues.xAT = fpeAdjustedValues.xAT;
kpeWithTendonDeltaValues.yAT = ...
    calcCentralDifferenceDataSeries(fpeAdjustedValues.xAT,...
                                    fpeAdjustedValues.yAT);


for i=1:1:length(fpe31Values.x)
    fiberKinematics = calcFixedWidthPennatedFiberKinematicsAlongTendon(...
                                                fpe31Values.x(i,1).*lceOpt,...
                                                0,...
                                                lceOpt,...
                                                alphaOpt);

    alpha = fiberKinematics.pennationAngle;
    lceAT = fiberKinematics.fiberLengthAlongTendon;


    fpeNAT = fiso*fpe31Values.y(i,1)*cos(alpha)/fisoAT;




    fpe31Values.xAT(i,1) = lceAT/lceOptAT;
    fpe31Values.yAT(i,1) = fpeNAT;

    if(fiberKinematics.isClamped)
        fprintf('%i. fal isClamped \n',i);
    end
end

disp('Checking the fpe and fpe31 produce similar values');
disp('Note: fpe31 is constructed using the ECM and titin curves');
for i=1:1:length(fpeValues.x)
    assert(abs(fpeValues.x(i,1)-fpe31Values.x(i,1)) < 1e-3,...
        ['Error: fpe fpe31 do not match for x at ',num2str(i)]);
    assert(abs(fpeValues.y(i,1)-fpe31Values.y(i,1)) < 5e-2,...
        ['Error: fpe fpe31 do not match for y at ',num2str(i)]);
    assert(abs(fpeValues.xAT(i,1)-fpe31Values.xAT(i,1)) < 1e-3,...
        ['Error: fpe fpe31 do not match for xAT at ',num2str(i)]);
    assert(abs(fpeValues.yAT(i,1)-fpe31Values.yAT(i,1)) < 5e-2,...
        ['Error: fpe fpe31 do not match for yAT at ',num2str(i)]);
end

%Evaluate the force-velocity curve along the tendon
for i=1:1:length(fvValues.x)
    fiberKinematics = calcFixedWidthPennatedFiberKinematicsAlongTendon(...
                                                lceOpt,...
                                                fvValues.x(i,1)*(vMax),...
                                                lceOpt,...
                                                alphaOpt);

    alpha = fiberKinematics.pennationAngle;
    lceAT = fiberKinematics.fiberLengthAlongTendon;
    vceAT = fiberKinematics.fiberVelocityAlongTendon;


    fvValues.xAT(i,1) = vceAT/vMaxAT;
    fvValues.yAT(i,1) = fiso*fvValues.y(i,1)*cos(alpha)/fisoAT;
end

fig=figure;
subplot(1,3,1);
    plot(falValues.x,falValues.y,'-','Color',[1,1,1].*0.5,'LineWidth',2);
    hold on;
    plot(falValues.xAT,falValues.yAT,'-','Color',[1,0,0]);
    hold on;
    
    xlabel('$$\tilde{\ell}^{M}$$');
    ylabel('$$\tilde{f}^{L}$$');
    box off;

subplot(1,3,2);
    plot(fpeValues.x,fpeValues.y,'--','Color',[1,1,1].*0.5,'LineWidth',2);
    hold on;
    plot(fpeValues.xAT,fpeValues.yAT,'-','Color',[0,0,1]);
    hold on;

    plot(fpe31Values.x,fpe31Values.y,'--','Color',[1,0.5,0.5]);
    hold on;    
    plot(fpe31Values.xAT,fpe31Values.yAT,'-','Color',[1,0,0]);
    hold on;



    plot(fpeAdjustedValues.x,fpeAdjustedValues.y,'--','Color',[0.5,0.5,1]);
    hold on;    
    plot(fpeAdjustedValues.xAT,fpeAdjustedValues.yAT,'-','Color',[0,0,1]);
    hold on;

    lceNfpe1 = interp1(fpeAdjustedValues.yAT,...
                       fpeAdjustedValues.xAT,...
                       1);
    kceNfpe1 = interp1(kpeWithTendonDeltaValues.xAT,...
                       kpeWithTendonDeltaValues.yAT,...
                       lceNfpe1);
    kcefpe1 = kceNfpe1*(1/1000)*(catSoleusHL1997.fceOptAT/catSoleusHL1997.lceOptAT);
    th=text(lceNfpe1,1,sprintf('%1.3f N/mm',kcefpe1),...
         'HorizontalAlignment','right',...
         'Color',[0,0,1]);
    th.Rotation=45;
    hold on;
   

    xlabel('$$\tilde{\ell}^{M}$$');
    ylabel('$$\tilde{f}^{PE}$$');
    box off;
    
subplot(1,3,3);
    plot(fvValues.x,fvValues.y,'-','Color',[0,0,0]);
    hold on;
    plot(fvValues.xAT,fvValues.yAT,'--','Color',[1,0,0]);
    hold on;
    
    xlabel('$$\tilde{v}^{M}$$');
    ylabel('$$\tilde{f}^{V}$$');
    box off;


fvNMax = calcQuadraticBezierYFcnXDerivative(1,...
       felineSoleusNormMuscleQuadraticCurves.fiberForceVelocityCurve,0);
fvNATMax    = fvNMax*cos(alphaOpt);


vNAT  = [fvValues.xAT; 1];
fvNAT = [fvValues.yAT; fvNATMax]; 

lceNAT_falN = falValues.xAT;
lceNAT_fpeN = fpeAdjustedValues.xAT;
vceNAT_fvNAT = vNAT;
dlceN1 = 0;

if(flag_solveForHerzogLeonard1997Params==1)
    dataHL1997=getHerzogLeonard1997Keypoints();

    fpeN0Exp = dataHL1997.fpe(1,1);
    lceN0Exp = dataHL1997.l(1,1);

    fpeN1Exp = dataHL1997.fpe(1,2);
    lceN1Exp = dataHL1997.l(1,2);

    idxMin = find(fpeAdjustedValues.yAT > 0.01,1);

    lceN1 = interp1(fpeAdjustedValues.yAT(idxMin:end,1),...
                    lceNAT_fpeN(idxMin:end,1),...
                    fpeN1Exp);
    dlceN1 = lceN1Exp-lceN1;

    fpeAdjustedValues.xAT=fpeAdjustedValues.xAT+dlceN1;

    kpeWithTendonDeltaValues.xAT = fpeAdjustedValues.xAT;
    kpeWithTendonDeltaValues.yAT = ...
        calcCentralDifferenceDataSeries(fpeAdjustedValues.xAT,...
                                        fpeAdjustedValues.yAT);
    
    subplot(1,3,2);
        plot(fpeAdjustedValues.xAT,...
             fpeAdjustedValues.yAT,'-c');
        hold on;
        plot(lceN1Exp,fpeN1Exp,'xk');
        hold on
        plot(lceN0Exp,fpeN0Exp,'xk');
        hold on

    
        lceNfpe1 = interp1(fpeAdjustedValues.yAT,...
                           fpeAdjustedValues.xAT,...
                           1);
        kceNfpe1 = interp1(kpeWithTendonDeltaValues.xAT,...
                           kpeWithTendonDeltaValues.yAT,...
                           lceNfpe1);
        kcefpe1 = kceNfpe1*(1/1000)*(catSoleusHL1997.fceOptAT/catSoleusHL1997.lceOptAT);
        th=text(lceNfpe1,1,sprintf('%1.3f N/mm',kcefpe1),...
             'HorizontalAlignment','right',...
             'Color',[0,1,1]);
        th.Rotation=45;

    hold on;        


   flN0Exp = dataHL1997.fa(1,1);
   flN1Exp = dataHL1997.fa(1,2);

   subplot(1,3,1);
        plot([0.5,1.5],[1,1].*flN0Exp,'-k');
        hold on
        plot(lceN0Exp,flN0Exp,'xk');
        hold on
        plot([0.5,1.5],[1,1].*flN1Exp,'-k');
        hold on     
        plot(lceN1Exp,flN1Exp,'xk');
        hold on
    
end

if(flag_adjustCurvesImplicitlyIncludeTendon==1)
    assert(0,'Error: this is not working at the moment');
    lceNAT_falN = (lceNAT_falN.*catSoleusHL2002.lceOpt) ...
        ./(catSoleusHL2002.lceOpt + catSoleusHL2002.ltSlk);

    lceNAT_fpeN = ((lceNAT_fpeN+dlceN1).*catSoleusHL2002.lceOpt) ...
        ./(catSoleusHL2002.lceOpt + catSoleusHL2002.ltSlk);

    vceNAT_fvNAT = (vceNAT_fvNAT.*(catSoleusHL2002.lceOpt*catSoleusHL2002.vceMax)) ...
        ./((catSoleusHL2002.lceOpt + catSoleusHL2002.ltSlk)*catSoleusHL2002.vceMax);

    subplot(1,3,1);
        plot(lceNAT_falN,falValues.yAT,'-','Color',[1,0,1]);
        hold on;
        
        xlabel('$$\tilde{\ell}^{M}$$');
        ylabel('$$\tilde{f}^{L}$$');
        box off;
    
    subplot(1,3,2);       
        plot(lceNAT_fpeN,fpeValues.yAT,'-','Color',[1,0,1]);
        hold on;
       
        xlabel('$$\tilde{\ell}^{M}$$');
        ylabel('$$\tilde{f}^{PE}$$');
        box off;
        
    subplot(1,3,3);
        plot(vceNAT_fvNAT,fvNAT,'--','Color',[1,0,1]);
        hold on;
        
        xlabel('$$\tilde{v}^{M}$$');
        ylabel('$$\tilde{f}^{V}$$');
        box off;

        
    success = writeFortranVector(lceNAT_falN, falValues.yAT, 10, ...
        'output/fortran/MAT156Tables/defaultFelineSoleusQ_activeForceLengthCurve_implicitRigidTendon.f');
    success = writeFortranVector(lceNAT_fpeN, fpeAdjustedValues.yAT, 11, ...
        'output/fortran/MAT156Tables/defaultFelineSoleusQ_forceLengthCurve_implicitRigidTendon.f');
    success = writeFortranVector(vceNAT_fvNAT, fvNAT, 12, ...
        'output/fortran/MAT156Tables/defaultFelineSoleusQ_forceVelocityCurve_implicitRigidTendon.f');
        
end

disp('To do: write umat43 and MAT 156 architectural properties to file');

success = writeFortranVector(falValues.xAT, falValues.yAT, 10, ...
    'output/fortran/MAT156Tables/defaultFelineSoleusQ_activeForceLengthCurve.f');
success = writeFortranVector(fpeAdjustedValues.xAT, fpeAdjustedValues.yAT, 11, ...
    'output/fortran/MAT156Tables/defaultFelineSoleusQ_forceLengthCurve.f');
success = writeFortranVector(vNAT, fvNAT, 12, ...
    'output/fortran/MAT156Tables/defaultFelineSoleusQ_forceVelocityCurve.f');


%%
% Remove the directories
%%
rmpath(parametersDirectoryTreeMTParams);
rmpath(parametersDirectoryTreeExperiments);
rmpath(parametersDirectoryTreeModels);
rmpath(parametersDirectoryTreeCurves);
