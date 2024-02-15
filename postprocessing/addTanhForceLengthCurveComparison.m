function figH = addTanhForceLengthCurveComparison(figH,...
    plotSettings, flag_usingOctave)

subPlotPanel            = plotSettings.subPlotPanel;
indexPlotRow            = plotSettings.indexPlotRow;
flag_plotBezierCurves   = plotSettings.flag_plotBezierCurves;
flag_plotTanhCurves     = plotSettings.flag_plotTanhCurves;
flag_plotTanCurves      = plotSettings.flag_plotTanCurves;
bezierColor             = plotSettings.bezierColor;
tanhColor               = plotSettings.tanhColor;
tanColor                = plotSettings.tanColor;

%%
% Create and evaluate the passive-force-length curves
%%
  flag_enableNumericallyNonZeroGradients  = 0;
  smallNumericallyNonZeroNumber           = sqrt(sqrt(eps));
  smallNumericallyNonZeroSlope            = sqrt(eps);

  normLengthZero = 1+0; 
  normLengthToe  = 1+0.6;
  fToe  = 1;
  yZero = 0;
  kZero = 0;
  if(flag_enableNumericallyNonZeroGradients)
    yZero   = smallNumericallyNonZeroValue;
    kZero   = smallNumericallyNonZeroSlope; 
  end         
  kLow  = 0.2;
  kToe  = 2/(normLengthToe-normLengthZero);
  curviness = 0.75;  
  flag_computeIntegral = 1;

  fiberForceLengthCurve = ...
    createFiberForceLengthCurve2021(normLengthZero,...
                                normLengthToe,...
                                fToe,...
                                yZero,...
                                kZero,...
                                kLow,...
                                kToe,...
                                curviness,...
                                flag_computeIntegral,...
                                'FiberForceLengthCurve',...
                                flag_usingOctave);   


x0      = fiberForceLengthCurve.xEnd(1,1);
x1      = fiberForceLengthCurve.xEnd(1,2);
y0      = fiberForceLengthCurve.yEnd(1,1);
y1      = fiberForceLengthCurve.yEnd(1,2);
dydx0   = fiberForceLengthCurve.dydxEnd(1,1);
dydx1   = fiberForceLengthCurve.dydxEnd(1,2);

yNegInf = y0;
yInf    = inf;  

xScale      = 0.9;
xAtIntYZero = 0;
xPoint = normLengthToe;
yPoint = fToe;

[A,B,C,D,E,F] = calcTanhSegmentCoefficientsUpd(...
                    x0,x1,dydx0,dydx1,...
                    yNegInf,yInf,...
                    xScale,xPoint, yPoint, xAtIntYZero);

%[A,B,C,D,E,F] = calcTanhSegmentCoefficients(x0,x1,dydx0,dydx1,...
%                                            yNegInf,yInf,...
%                                            xShift,xScale,xAtIntYZero);
forceLengthTanhCoeffs = [A,B,C,D,E,F];  

npts = 100;
lceN = [0:(normLengthToe+0.1)/(npts-1):(normLengthToe+0.1)]';

fpeBezierSample = zeros(npts,3);
fpeTanhSample   = zeros(npts,3);

for i=1:1:3
    for j=1:1:npts
        fpeBezierSample(j,i) = calcBezierYFcnXDerivative(lceN(j,1),...
                                              fiberForceLengthCurve,i-2);

        fpeTanhSample(j,i) = calcTanhSeriesDerivative(lceN(j,1),...
                                       forceLengthTanhCoeffs,i-2);
    end

    figure(figH);
    subplot('Position',reshape(subPlotPanel(2,i,:),1,4));

    if(flag_plotBezierCurves==1)
        plot( lceN,fpeBezierSample(:,i),...
              'Color',bezierColor,'LineWidth',1,...
              'DisplayName','Bezier');
        hold on;
    end
    if(flag_plotTanhCurves==1)
        plot( lceN,fpeTanhSample(:,i),...
              'Color',tanhColor,'LineWidth',1,...
              'DisplayName','Tanh');
        hold on;
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
            title('Force-Length Curve Integral');
            xlabel('Norm. Length ($\ell/\ell^M_o$)');
            ylabel('Norm. Energy ($\tilde{f}^M \, \tilde{\ell}^M$)');
            legend('Location','NorthEast');            
            legend boxoff;
        case 2
            plot([x0;x1],[y0;y1],'o','MarkerSize',5,...
                'Color',bezierColor,'MarkerFaceColor',[1,1,1]);
            hold on;
            title('Force-Length  Curve Value');            
            xlabel('Norm. Length ($\ell/\ell^M_o$)')
            ylabel('Norm. Force ($f/f^M_o$)');            

            
        case 3
            plot([x0;x1],[dydx0;dydx1],'o','MarkerSize',5,...
                'Color',bezierColor,'MarkerFaceColor',[1,1,1]);
            hold on;

            title('Derivative');                        
            xlabel('Norm. Length ($\ell/\ell^M_o$)')
            ylabel('Norm. Stiffness ($\delta \tilde{f^M} / \delta \tilde{\ell^M}$)');
            hold on;
            
        otherwise
            assert(0,'Error: missing postprocessing code for the current derivative');
    end

end
