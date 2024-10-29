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
catSoleusHL2002.vceMax  =4.650;

catSoleusHL1997.lceOpt  =0.038;
catSoleusHL1997.fceOpt  =41.661837;
catSoleusHL1997.lceOptAT=0.0377168;
catSoleusHL1997.fceOptAT=41.351296;
catSoleusHL1997.lmtOptAT=0.0681679;
catSoleusHL1997.penOpt  =0.1221730;
catSoleusHL1997.penOptD =7.0 ;
catSoleusHL1997.ltSlk   =0.0304511;
catSoleusHL1997.et      =0.0458333;
catSoleusHL1997.vceMax  =4.650;


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

lceN2Exp = 0.792;
flN2Exp  = 0.88-(0.897-0.91);


lceOpt  = catSoleusHL1997.lceOptAT;
dWdes   = 0.470961;
nuCEdes = 1.328409;
Fmax    = 41.351296;
FPEE    = 2.031664339235247;
FPEEupd = FPEE*1.;
LPEE0   = 0.817625283209556;
LPEE0upd =LPEE0; 
nuPEE   = 1.5907747; 


%%
%Sample tendon force length curve
%%
lSEE0    = 0.0304511;
dUSEEnll = 0.0347222;
dUSEEl   = 0.0222222;
dFSEE0   = 27.5735;%30*(30/32.64);%14.337374;

ltN = [1:0.001:1.12]';
fsee = zeros(length(ltN),1);
for i=1:1:length(ltN)
    fsee(i,1)=calcFseeUmat41(lSEE0*ltN(i,1),lSEE0,dUSEEnll,dUSEEl,dFSEE0);
end
fseeN=fsee./Fmax;
idxMin = find(fseeN > 0.01);
ltNFtFiso = interp1(fseeN(idxMin:end,1),ltN(idxMin:end,1),1);

%idxFtFmax =  find(fseeN>1,1);
idxLtMin = find(fseeN > 0.01,1);
ltNFmax = interp1(fseeN(idxLtMin:end),ltN(idxLtMin:end),1);
kseN = calcCentralDifferenceDataSeries(ltN,fseeN);
kseNFiso = interp1(ltN(idxLtMin:end),kseN(idxLtMin:end),ltNFmax);

%%
% Solve for the parameters that shift the fpe curve to the left until
% it intersects the desired point
%%

fpeeParams.lceOpt   = lceOpt;
fpeeParams.dWdes    = dWdes;
fpeeParams.nuCEdes  = nuCEdes;
fpeeParams.Fmax     = Fmax;
fpeeParams.FPEE     = FPEE;
fpeeParams.LPEE0    = LPEE0;
fpeeParams.nuPEE    = nuPEE;
fpeeParams.lSEE0    = lSEE0;
fpeeParams.dUSEEnll = dUSEEnll;
fpeeParams.dUSEEl   = dUSEEl;
fpeeParams.dFSEE0   = dFSEE0;

lceFmax = calcFpeeInverseUmat41(Fmax, lceOpt,dWdes,Fmax,FPEE,LPEE0,nuPEE);
dfpee   = calcFpeeDerivativeUmat41(lceFmax, lceOpt,dWdes,Fmax,FPEE,LPEE0,nuPEE);

mm2m=0.001;
manualScalingKmt = 1;%(3.78/3.88);

% Target of 3.7814 N/mm comes from wanting to match the stiffness of
% the umat43 (VEXAT) cat soleus model
%kmtTarget = ;
%kpeTarget = ((1/kmtTarget)-(1/kseNFiso))^-1;

%ltN1Exp = calcFseeInverseUmat41(fpeN1Exp,lSEE0,dUSEEnll,dUSEEl,dFSEE0);
%dlceN1ltN1 = (ltN1Exp-lSEE0)/lceOpt;


targetPointsFpe.lceN0  = lceN1Exp;
targetPointsFpe.fpeN0  = fpeN1Exp;
targetPointsFpe.fpeN1  = Fmax/Fmax;
targetPointsFpe.dfpeN1 = 3.7814*manualScalingKmt*(1/mm2m)*(lceOpt/Fmax);


errFcnFpe = @(arg)calcFpeErrEHTMM(arg, fpeeParams,targetPointsFpe);
x0 = [LPEE0upd;FPEEupd];
[x, fval] = lsqnonlin(errFcnFpe,x0);
LPEE0upd=x(1,1);
FPEEupd=x(2,1);

disp('Optimized fpe parameters');
fprintf('%e\tLPEE0\n%e\tFPEE\n',LPEE0upd,FPEEupd);

%%
% Active force-length curve
%%

dWdes   = 0.470961;
nuCEdes = 1.328409;
dWasc   = 0.5656740;%0.4005710;
nuCEasc = 2.082109;%4.5843765;

fisomParams.dWdes   = dWdes;
fisomParams.nuCEdes = nuCEdes;
fisomParams.dWasc   = dWasc;
fisomParams.nuCEasc = nuCEasc;
fisomParams.lceOpt=lceOpt;

targetPointsFisom.lceN0   = lceN0Exp;
targetPointsFisom.flN0    = flN0Exp;
targetPointsFisom.lceN1   = lceN1Exp;
targetPointsFisom.flN1    = flN1Exp;

targetPointsFisom.lceN2   = lceN2Exp;
targetPointsFisom.flN2    = flN2Exp;


errFcnFisom= @(arg)calcFisomErrEHTMM(arg, fisomParams,targetPointsFisom);
x0 = [dWasc;nuCEasc];
[x, fval] = lsqnonlin(errFcnFisom,x0);
dWasc   = x(1,1);
nuCEasc = x(2,1);

disp('Optimized fisom parameters');
fprintf('%e\tdWasc\n%e\tnuCEasc\n',dWasc,nuCEasc);

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
    subplot(1,4,1);
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
        plot([0.5,1.5],[1,1].*flN2Exp,'-k');
        hold on     
        plot(lceN2Exp,flN2Exp,'xk');
        hold on

        xlabel('Norm. Length ($$\ell/\ell^M_o$$)');
        ylabel('Norm. Force ($$f/f^M_o$$)');
        title('Active Force Length Relation');
    
    subplot(1,4,2);
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

    subplot(1,4,3);

        dfpeNum = calcCentralDifferenceDataSeries(lceN,fpeN);

        plot(lceN,dfpeN,'-k');
        hold on;
        plot(lceN,dfpeNum,'--r');

        xlabel('Norm. Length ($$\ell/\ell^M_o$$)');
        ylabel('Norm. Stiffness ($$(f/f^M_o)/\ell$$)');
        title('Passive Stiffness Length Relation');

    subplot(1,4,4);
        plot(ltN,fseeN);
        hold on;
        text(ltNFtFiso,1,sprintf('eT=%1.6f\nkT=%1.6f',ltNFtFiso,kseNFiso));

        xlabel('Norm. Length ($$\ell^T/\ell^T_{s}$$)');
        ylabel('Norm. Force ($$f^T/f^M_o$$)');
        title('Tendon Force Length');