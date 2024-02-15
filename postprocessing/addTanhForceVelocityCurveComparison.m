function figH = addTanhForceVelocityCurveComparison(figH,...
    plotSettings, flag_usingOctave)

flag_plotComponents=0;
subPlotPanel            = plotSettings.subPlotPanel;
indexPlotRow            = plotSettings.indexPlotRow;
flag_plotBezierCurves   = plotSettings.flag_plotBezierCurves;
flag_plotTanhCurves     = plotSettings.flag_plotTanhCurves;
flag_plotTanCurves      = plotSettings.flag_plotTanCurves;
bezierColor             = plotSettings.bezierColor;
tanhColor               = plotSettings.tanhColor;
tanColor                = plotSettings.tanColor;


flag_enableNumericallyNonZeroGradients  = 0;
smallNumericallyNonZeroNumber           = sqrt(sqrt(eps));
smallNumericallyNonZeroSlope            = sqrt(eps);

curvinessEccentricForceVelocity = 1.0;
flag_sharpEccentricTransition = 0;

forceVelocityMultiplierAtHalfMaximumFiberVelocity       = 0.2;
forceVelocityMultiplierAtLowEccentricFiberVelocity      = 1.4;
forceVelocityMultiplierAtMaximumEccentricFiberVelocity  = 1.6;

fiberForceVelocityCurve ...
  = createFiberForceVelocityCurve(...
      forceVelocityMultiplierAtHalfMaximumFiberVelocity,...
      forceVelocityMultiplierAtLowEccentricFiberVelocity,...
      forceVelocityMultiplierAtMaximumEccentricFiberVelocity,...
      curvinessEccentricForceVelocity,...
      flag_sharpEccentricTransition,...
      flag_enableNumericallyNonZeroGradients,...
      smallNumericallyNonZeroNumber,...
      smallNumericallyNonZeroSlope,...
      'ForceVelocityCurve',...
      flag_usingOctave);



x0   = -1.05;
y0   = calcBezierYFcnXDerivative(x0,fiberForceVelocityCurve,0);
dydx0= 0;%calcBezierYFcnXDerivative(x0,fiberForceVelocityCurve,1);

x1   = 0;
y1   = calcBezierYFcnXDerivative(x1,fiberForceVelocityCurve,0);
dydx1= calcBezierYFcnXDerivative(x1,fiberForceVelocityCurve,1);


yNegInf = 0;
yInf    = inf;  

xScale      = 1.75;%1.75;
xAtIntYZero = x0;
xPoint = x1;
yPoint = y1;

[A,B,C,D,E,F] = calcTanhSegmentCoefficientsUpd(...
                    x0,x1,dydx0,dydx1,...
                    yNegInf,yInf,...
                    xScale,xPoint, yPoint, xAtIntYZero);

%[A,B,C,D,E,F] = calcTanhSegmentCoefficients(x0,x1,dydx0,dydx1,...
%                                            yNegInf,yInf,...
%                                            xShift,xScale,xAtIntYZero);
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

%Add a curve that will flatten the eccentric side to the desired final
%slope
dydx1 = (forceVelocityMultiplierAtMaximumEccentricFiberVelocity...
       -forceVelocityMultiplierAtLowEccentricFiberVelocity)/1 -dy;
xAtIntYZero=0;
xScale = 1;
[A,B,C,D,E,F] = calcTanhSegmentCoefficientsUpd(...
                    0.0,0.08,0,dydx1,...
                    yNegInf,yInf,...
                    xScale,0, 0, xAtIntYZero);
forceVelocityTanhCoeffs = [forceVelocityTanhCoeffs;... 
                           A,B,C,D,E,F];




npts = 100;
xP0 = -1.1;
xP1 = 1.3;
vceN = [(xP0):(xP1-xP0)/(npts-1):(xP1)]';

fvBezierSample = zeros(npts,3);

fvTanhSample        = zeros(npts,3);
fvComponentSample   = zeros(npts,size(forceVelocityTanhCoeffs,1));

for i=2:1:3
    for j=1:1:npts
        fvBezierSample(j,i) = calcBezierYFcnXDerivative(vceN(j,1),...
                                              fiberForceVelocityCurve,i-2);

        fvTanhSample(j,i) = calcTanhSeriesDerivative(vceN(j,1),...
                                       forceVelocityTanhCoeffs,i-2);

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
%     if(flag_plotTanCurves==1)
%         plot( lceN,fpeTanSample(:,i),...
%               'Color',tanColor,'LineWidth',1,...
%               'DisplayName','Tan');
%         hold on;
%     end
    box off;

    switch(i)
        case 1
            title('Force-Velocity Curve Integral');
            xlabel('Norm. Velocity ($v/v^M_{max}$)');
            ylabel('Norm. Int. Force-Velocity ($\int \tilde{f}^M \tilde{v}^M$)');
            legend('Location','NorthEast');            
            legend boxoff;
        case 2
            plot([x0;x1],[y0;y1],'o','MarkerSize',5,...
                'Color',bezierColor,'MarkerFaceColor',[1,1,1]);
            hold on;
            title('Force-Velocity Curve Value');            
            xlabel('Norm. Velocity ($v/v^M_{max}$)')
            ylabel('Norm. Force ($f/f^M_o$)');            

            
        case 3
            plot([x0;x1],[dydx0;dydx1],'o','MarkerSize',5,...
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