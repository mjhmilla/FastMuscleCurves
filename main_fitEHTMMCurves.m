clc;
close all;
clear all;


flag_solveForHerzogLeonard1997Params=1;

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
%
%%

dataHL1997=getHerzogLeonard1997Keypoints();
fpeN0Exp = dataHL1997.fpe(1,1);
lceN0Exp = dataHL1997.l(1,1);

fpeN1Exp = dataHL1997.fpe(1,2);
lceN1Exp = dataHL1997.l(1,2);


flN0Exp = dataHL1997.fa(1,1);
flN1Exp = dataHL1997.fa(1,2);

lceOpt  = catSoleusHL1997.lceOptAT;
dWdes   = 0.3852852;
nuCEdes = 1.9043187;
Fmax    = 41.351296;
FPEE    = 0.654342;
FPEEupd = FPEE*1.7;
LPEE0   = 1.0063986;
LPEE0upd =1.0063986-0.203145; 
nuPEE   = 1.5907747; 

%%
% Solve for the parameters that shift the fpe curve to the left until
% it intersects the desired point
%%

fpeeParams.lceOpt=lceOpt;
fpeeParams.dWdes=dWdes;
fpeeParams.nuCEdes=nuCEdes;
fpeeParams.Fmax=Fmax;
fpeeParams.FPEE=FPEE;
fpeeParams.LPEE0=LPEE0;
fpeeParams.nuPEE=nuPEE;

lceFmax = calcFpeeInverseUmat41(Fmax, lceOpt,dWdes,Fmax,FPEE,LPEE0,nuPEE);
dfpee   = calcFpeeDerivativeUmat41(lceFmax, lceOpt,dWdes,Fmax,FPEE,LPEE0,nuPEE);

targetPoints.lceN0 = lceN1Exp;
targetPoints.fpeN0 = fpeN1Exp;
targetPoints.lceN1 = lceFmax/lceOpt;
targetPoints.fpeN1 = Fmax/Fmax;
targetPoints.dfpeN1 = dfpee*(lceOpt/Fmax);


errFcnFpe = @(arg)calcFpeErrEHTMM(arg, fpeeParams,targetPoints);
x0 = [LPEE0upd;FPEEupd];
[x, fval] = fminsearch(errFcnFpe,x0);
LPEE0upd=x(1,1);
FPEEupd=x(2,1);

disp('Optimized fpe parameters');
fprintf('%e\tLPEE0\n%e\tFPEE',LPEE0upd,FPEEupd);

%%
% Active force-length curve
%%

dWdes   = 0.3852852;
nuCEdes = 1.9043187;
dWasc   = 0.57;%0.4005710;
nuCEasc = 1.9;%4.5843765;

fisomParams.dWdes   = dWdes;
fisomParams.nuCEdes = nuCEdes;
fisomParams.dWasc   = dWasc;
fisomParams.nuCEasc = nuCEasc;
fisomParams.lceOpt=lceOpt;

targetPoints.lceN0   = lceN0Exp;
targetPoints.flN0    = flN0Exp;
targetPoints.lceN1   = lceN1Exp;
targetPoints.flN1    = flN1Exp;

errFcnFisom= @(arg)calcFisomErrEHTMM(arg, fisomParams,targetPoints);
x0 = [dWasc;nuCEasc];
[x, fval] = fminsearch(errFcnFisom,x0);
dWasc   = x(1,1);
nuCEasc = x(2,1);

disp('Optimized fisom parameters');
fprintf('%e\tdWasc\n%e\tnuCEasc',dWasc,nuCEasc);

%%
% Solve for the dWasc and nuCEasc parameters that intersect the desired 
% points
%%

lceN = [0:0.01:1.6]';
fpeN = zeros(size(lceN));
dfpeN = zeros(size(lceN));

fpeNupd = zeros(size(lceN));

falN = zeros(size(lceN));

for i=1:1:length(lceN)
    fpeN(i,1) = calcFpeeUmat41(lceN(i,1)*lceOpt, lceOpt,dWdes,...
                               Fmax,FPEE,LPEE0,nuPEE);
    fpeN(i,1) = fpeN(i,1)/Fmax;

    dfpeN(i,1)= calcFpeeDerivativeUmat41(lceN(i,1)*lceOpt, lceOpt,dWdes,Fmax,FPEE,LPEE0,nuPEE);
    dfpeN(i,1)=dfpeN(i,1)*(1/Fmax)*(lceOpt);

    fpeNupd(i,1) = calcFpeeUmat41(lceN(i,1)*lceOpt, lceOpt,dWdes,...
                               Fmax,FPEEupd,LPEE0upd,nuPEE);
    fpeNupd(i,1) = fpeNupd(i,1)/Fmax;

    falN(i,1) = calcFisomUmat41(lceN(i,1)*lceOpt,lceOpt,dWdes, nuCEdes, dWasc, nuCEasc);
end


idxMin = find(fpeN > 0.01,1);

lceN1 = interp1(fpeN(idxMin:end,1),...
                lceN(idxMin:end,1),...
                fpeN1Exp);
dlceN1 = lceN1Exp-lceN1;

idxMin = find(fpeNupd > 0.01,1);

lceN1upd = interp1(fpeNupd(idxMin:end,1),...
                lceN(idxMin:end,1),...
                fpeN1Exp);
dlceN1upd = lceN1Exp-lceN1upd;


fig=figure;
    subplot(1,3,1);
        plot(lceN,falN,'-k');
        hold on;
        plot([0.5,1.5],[1,1].*flN0Exp,'-k');
        hold on
        plot(lceN0Exp,flN0Exp,'xk');
        hold on
        plot([0.5,1.5],[1,1].*flN1Exp,'-k');
        hold on     
        plot(lceN1Exp,flN1Exp,'xk');
        hold on

        xlabel('Norm. Length ($$\ell/\ell^M_o$$)');
        ylabel('Norm. Force ($$f/f^M_o$$)');
        title('Active Force Length Relation');
    
    subplot(1,3,2);
        plot(lceN,fpeN,'-k');
        hold on;
        plot(lceN,fpeNupd,'-b');
        hold on;
        plot(lceN1Exp,fpeN1Exp,'xk');
        hold on
        plot(lceN0Exp,fpeN0Exp,'xk');
        hold on
        text(lceN1Exp,fpeN1Exp+0.05,...
            sprintf('%1.6f',lceN1Exp),...
            'HorizontalAlignment','right');
        hold on
        text(lceN1,fpeN1Exp+0.05,...
            sprintf('%1.6f',lceN1),...
            'HorizontalAlignment','left');  
        text(lceN1,fpeN1Exp-0.05,...
            sprintf('%1.6f',dlceN1),...
            'HorizontalAlignment','center');  
        hold on;
        text(lceN1upd,fpeN1Exp-0.05,...
            sprintf('%1.6f',dlceN1upd),...
            'HorizontalAlignment','center');  
        hold on;

        plot([lceN1Exp,lceN1],...
            [1,1].*fpeN1Exp,'-r');

        xlabel('Norm. Length ($$\ell/\ell^M_o$$)');
        ylabel('Norm. Force ($$f/f^M_o$$)');
        title('Passive Force Length Relation');

    subplot(1,3,3);

        dfpeNum = calcCentralDifferenceDataSeries(lceN,fpeN);

        plot(lceN,dfpeN,'-k');
        hold on;
        plot(lceN,dfpeNum,'--r');

        xlabel('Norm. Length ($$\ell/\ell^M_o$$)');
        ylabel('Norm. Stiffness ($$(f/f^M_o)/\ell$$)');
        title('Passive Stiffness Length Relation');
