function figH = addTanhForceVelocityCurveComparison(figH,...
    fiberForceVelocityCurve,...
    fvDomainTest,...
    plotSettings)

flag_fitfvCurve=1;
flag_plotComponents=0;
subPlotPanel            = plotSettings.subPlotPanel;
indexPlotRow            = plotSettings.indexPlotRow;
flag_plotBezierCurves   = plotSettings.flag_plotBezierCurves;
flag_plotTanhCurves     = plotSettings.flag_plotTanhCurves;
flag_plotTanCurves      = plotSettings.flag_plotTanCurves;
bezierColor             = plotSettings.bezierColor;
tanhColor               = plotSettings.tanhColor;
tanhErrorColor          = plotSettings.tanhErrorColor;
tanColor                = plotSettings.tanColor;



tanhSeriesParams(2) = ...
    struct('x0',0,'dydx0',0,...
           'x1',0,'dydx1',0,...
           'yNegInf',0,'yInf',0,...
           'xScale',0,'xPoint',0,...
           'yPoint',0,'xAtIntYZero',0);

optParams(2) = struct('names',{''});

yNegInf = 0;
yInf    = inf;  
dyPoint = 0.01;

tanhSeriesParams(1).x0          = -2.5;
tanhSeriesParams(1).x1          = 0;
tanhSeriesParams(1).dydx0       = 0;
tanhSeriesParams(1).dydx1       = calcBezierYFcnXDerivative(...
                                    tanhSeriesParams(1).x1,...
                                    fiberForceVelocityCurve,1)*1.1;
tanhSeriesParams(1).yNegInf     = yNegInf;
tanhSeriesParams(1).yInf        = yInf;
tanhSeriesParams(1).xScale      = 1.0;
tanhSeriesParams(1).xPoint      = 0;
tanhSeriesParams(1).yPoint      = 1+dyPoint;
tanhSeriesParams(1).xAtIntYZero = tanhSeriesParams(1).x0;

optParams(1).names = {'x0','dydx1','xScale'};
args = [tanhSeriesParams(1).x0;...
        tanhSeriesParams(1).dydx1;...
        tanhSeriesParams(1).xScale];

[A,B,C,D,E,F] = calcTanhSegmentCoefficientsUpd( ...
                    tanhSeriesParams(1).x0,     ...
                    tanhSeriesParams(1).x1,     ...
                    tanhSeriesParams(1).dydx0,  ...
                    tanhSeriesParams(1).dydx1,  ...
                    tanhSeriesParams(1).yNegInf,...
                    tanhSeriesParams(1).yInf,   ...
                    tanhSeriesParams(1).xScale, ...
                    tanhSeriesParams(1).xPoint, ...
                    tanhSeriesParams(1).yPoint, ...
                    tanhSeriesParams(1).xAtIntYZero);

forceVelocityTanhCoeffs = [A,B,C,D,E,F];  

%Get the slope at postive infinity
dy=0;
for i=1:1:size(forceVelocityTanhCoeffs,1)
    A = forceVelocityTanhCoeffs(i,1);
    B = forceVelocityTanhCoeffs(i,2);
    C = forceVelocityTanhCoeffs(i,3);
    D = forceVelocityTanhCoeffs(i,4);            
    dy = dy + (A + D);
end



tanhSeriesParams(2).x0          =0;
tanhSeriesParams(2).x1          =0.5;
tanhSeriesParams(2).dydx0       =0;
tanhSeriesParams(2).dydx1       = ...
    calcBezierYFcnXDerivative(1,fiberForceVelocityCurve,1) ...
   -calcTanhSeriesDerivative(1e6,forceVelocityTanhCoeffs,1);
tanhSeriesParams(2).yNegInf     =yNegInf;
tanhSeriesParams(2).yInf        =sign(tanhSeriesParams(2).dydx1)*yInf;
tanhSeriesParams(2).xScale      =1;
tanhSeriesParams(2).xPoint      =0;
tanhSeriesParams(2).yPoint      =-dyPoint;
tanhSeriesParams(2).xAtIntYZero =0;

optParams(2).names = {'x1','dydx1','xScale'};

args        = [args; ...
              tanhSeriesParams(2).x1; ...
              tanhSeriesParams(2).dydx1; ...
              tanhSeriesParams(2).xScale];
argsScaling = args;
argsScaled  = ones(size(args));


[A,B,C,D,E,F] = calcTanhSegmentCoefficientsUpd( ...
                    tanhSeriesParams(2).x0,     ...
                    tanhSeriesParams(2).x1,     ...
                    tanhSeriesParams(2).dydx0,  ...
                    tanhSeriesParams(2).dydx1,  ...
                    tanhSeriesParams(2).yNegInf,...
                    tanhSeriesParams(2).yInf,   ...
                    tanhSeriesParams(2).xScale, ...
                    tanhSeriesParams(2).xPoint, ...
                    tanhSeriesParams(2).yPoint, ...
                    tanhSeriesParams(2).xAtIntYZero);

forceVelocityTanhCoeffs = [forceVelocityTanhCoeffs;... 
                           A,B,C,D,E,F];

xA = fiberForceVelocityCurve.xEnd(1,1);
xB = 0;

nptsConc = 50;
concDomain = [xA:((xB-xA)/(nptsConc-1)):xB]';

xC = 0.01;
xD = 1;
nptsEccSlow = 20;
eccDomainSlow = [xC:((xD-xC)/(nptsEccSlow-1)):xD]';

xE = fiberForceVelocityCurve.xEnd(1,2)*2;
xF = fiberForceVelocityCurve.xEnd(1,2)*3;
nptsEccFast = 5;
eccDomainFast = [xE:((xF-xE)/(nptsEccFast-1)):xF]';


optDomain = [concDomain; ...
             eccDomainSlow; ...
             eccDomainFast];


errVec0 = calcTanhCurveError(argsScaled, optParams,...
            tanhSeriesParams,fiberForceVelocityCurve, ...
            optDomain,argsScaling);



if(flag_fitfvCurve==1)

    errFcn = @(argInput)calcTanhCurveError(argInput,...
                 optParams,tanhSeriesParams,fiberForceVelocityCurve,...
                 optDomain,argsScaling);

    [argScaledUpd,resnorm,residual,exitflag,output]=...
        lsqnonlin(errFcn,argsScaled);

    argUpd = argScaledUpd.*argsScaling;
%     argsScaling = argUpd;
%     argUpd = ones(size(argUpd));
% 
%     errFcn = @(argInput)calcTanhCurveError(argInput,...
%                  optParams,tanhSeriesParams,fiberForceVelocityCurve,...
%                  optDomain,argsScaling);
% 
%     [argScaledUpd,resnorm,residual,exitflag,output]=...
%         lsqnonlin(errFcn,argUpd);
%     
%     argUpd = argScaledUpd.*argsScaling;

    errVec1 = calcTanhCurveError(argScaledUpd, optParams,...
                tanhSeriesParams,fiberForceVelocityCurve,...
                optDomain,argsScaling);
    




    fprintf('%1.2e\tStarting Error\n%1.2e\tEnding Error\n',...
             sqrt(sum(errVec0.^2)),sqrt(sum(errVec1.^2)));  
    fprintf('%i\tExit flag',exitflag);

    localParams=tanhSeriesParams;
    idx=1;
    for i=1:1:length(optParams)
    
        varNames = optParams(i).names;
        for j=1:1:length(varNames)
            localParams(i).(varNames{j})=argUpd(idx,1);
            idx=idx+1;
        end
    
        x0_          =localParams(i).x0;
        x1_          =localParams(i).x1;
        dydx0_       =localParams(i).dydx0;
        dydx1_       =localParams(i).dydx1;
        yNegInf_     =localParams(i).yNegInf;
        yInf_        =localParams(i).yInf;
        xScale_      =localParams(i).xScale;
        xPoint_      =localParams(i).xPoint;
        yPoint_      =localParams(i).yPoint;
        xAtIntYZero_ =localParams(i).xAtIntYZero;
    
        [A,B,C,D,E,F] = calcTanhSegmentCoefficientsUpd(...
                        x0_,x1_,dydx0_,dydx1_,...
                        yNegInf_,yInf_,...
                        xScale_,xPoint_, yPoint_, xAtIntYZero_);
    
        forceVelocityTanhCoeffs(i,:) = [A,B,C,D,E,F];
    end
end

here=1;


vceN = fvDomainTest;
npts = length(fvDomainTest);

idxCE = find(vceN >= -1 & vceN <= 1);
idxC = find(vceN >= -1 & vceN <= 0);
idxE = find(vceN >= 0 & vceN <= 1);

fvBezierSample = zeros(npts,3);
fvErrorSample = zeros(npts,3);

fvTanhSample        = zeros(npts,3);
fvComponentSample   = zeros(npts,size(forceVelocityTanhCoeffs,1));

for i=2:1:3
    for j=1:1:npts
        fvBezierSample(j,i) = calcBezierYFcnXDerivative(vceN(j,1),...
                                              fiberForceVelocityCurve,i-2);

        fvTanhSample(j,i) = calcTanhSeriesDerivative(vceN(j,1),...
                                       forceVelocityTanhCoeffs,i-2);

%         if(i==2)
%             fvTanhSample(j,i)=fvTanhSample(j,i) + 0.1 + 0.1*tanh( (vceN(j,1)-0.03)/0.01 );
%         end

        fvErrorSample(j,i)=fvTanhSample(j,i)-fvBezierSample(j,i);        

%         if(i==2)
%             y0 =0.1+0.1*tanh((0-0.05)/0.01);
%             fvTanhSample(j,i) = fvTanhSample(j,i)+ 0.1+0.1*tanh((vceN(j,1)-0.05)/0.01)-y0;
%         end

        if(i-2 ==0)
            for k=1:1:size(fvComponentSample,2)
                fvComponentSample(j,k) = calcTanhSeriesDerivative(vceN(j,1),...
                                       forceVelocityTanhCoeffs(k,:),i-2);
            end
        end
    end

    figure(figH);
    subplot('Position',reshape(subPlotPanel(indexPlotRow ,i,:),1,4));

    if(flag_plotBezierCurves==1)
        plot( vceN,fvBezierSample(:,i),...
              'Color',bezierColor,'LineWidth',1,...
              'DisplayName','Bezier');
        hold on;
    end
    if(flag_plotTanhCurves==1)
        plot( vceN,fvTanhSample(:,i),...
              'Color',tanhColor,'LineWidth',1,...
              'DisplayName','Tanh');
        hold on;
        if(i==2 && flag_plotComponents==1)
            for k=1:1:size(fvComponentSample,2)
                plot( vceN,fvComponentSample(:,k),...
                      'Color',tanhColor.*0.5+[1,1,1].*0.5,'LineWidth',1,...
                      'DisplayName','Tanh-Component');
                hold on;
            end
        end

    end
    if(flag_plotBezierCurves==1 && flag_plotTanhCurves==1)
        plot(vceN,fvErrorSample(:,i),'-',...
            'Color',tanhErrorColor,...
            'DisplayName','Tanh Error');
        hold on;

        text(-1,1,sprintf('%1.2e: RMSE\n%1.2e: Conc. RMSE\n%1.2e: Ecc. RMSE\n',...
                  sqrt(sum(fvErrorSample(idxCE,i).^2)),...
                  sqrt(sum(fvErrorSample(idxC,i).^2)),...
                  sqrt(sum(fvErrorSample(idxE,i).^2)) ),...
                  'FontSize',6,...
                  'Color',tanhColor);
        hold on;

    end

    box off;

    switch(i)
        case 1
            title('Force-Velocity Curve Integral');
            xlabel('Norm. Velocity ($v/v^M_{max}$)');
            ylabel('Norm. Int. Force-Velocity ($\int \tilde{f}^M \tilde{v}^M$)');
            legend('Location','NorthEast');            
            legend boxoff;
        case 2

            y0=calcTanhSeriesDerivative(tanhSeriesParams(1).x0,...
                                       forceVelocityTanhCoeffs,i-2);
            y1=calcTanhSeriesDerivative(tanhSeriesParams(1).x1,...
                                       forceVelocityTanhCoeffs,i-2);
            
            plot([tanhSeriesParams(1).x0;...
                  tanhSeriesParams(1).x1],...
                  [y0;y1],'o','MarkerSize',5,...
                'Color',bezierColor,'MarkerFaceColor',[1,1,1]);
            hold on;
            title('Force-Velocity Curve Value');            
            xlabel('Norm. Velocity ($v/v^M_{max}$)')
            ylabel('Norm. Force ($f/f^M_o$)');            

            
        case 3

            dydx0=calcTanhSeriesDerivative(tanhSeriesParams(1).x0,...
                                       forceVelocityTanhCoeffs,i-2);
            dydx1=calcTanhSeriesDerivative(tanhSeriesParams(1).x1,...
                                       forceVelocityTanhCoeffs,i-2);            
            plot([tanhSeriesParams(1).x0;...
                  tanhSeriesParams(1).x1],...
                  [dydx0;dydx1],'o','MarkerSize',5,...
                'Color',bezierColor,'MarkerFaceColor',[1,1,1]);
            hold on;

            title('Derivative');                        
            xlabel('Norm. Velocity ($v/v^M_{max}$)')
            ylabel('Norm. Slope ($(f/f^M_o)/(v/v^M_{max})$)');
            hold on;
            
        otherwise
            assert(0,'Error: missing postprocessing code for the current derivative');
    end

end