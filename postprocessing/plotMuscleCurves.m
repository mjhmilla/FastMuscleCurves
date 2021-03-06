function [fig_pubCurves,subPlotPanel] = plotMuscleCurves(...
                      felineSoleusNormMuscleCurves,...
                      activeForceLengthCurveAnnotationPoints,...
                      felineSoleusMusculotendonPropertiesUpd,...
                      felineSoleusSarcomerePropertiesUpd,...
                      felineSoleusActiveForceLengthDataUpd,...
                      felineSoleusPassiveForceLengthDataUpd,...
                      pubOutputFolder)
                                    
  fig_pubCurves=figure;
  
  numberOfHorizontalPlotColumns = 3;
  numberOfVerticalPlotRows      = 2;

  pageWidth         = 21;
  pageHeight        = 29.7;
  plotWidth         = 4.5;
  plotHeight        = 4.5;
  plotHorizMarginCm = 1.5;
  plotVertMarginCm  = 2.;

  flag_usingOctave  = 0;
  subPlotPanel = [];
  plotConfigGeneric;  
    
  lineColorTendon = [0,0,0];
  lineWidthTendon = 1;
  
  labelRotationOffset = -10;
  
  plotProps(6) = struct('xlim',[],'ylim',[],'domain',[],...
                     'xticks',[],'xticklabels',{''},...
                     'yticks',[],'yticklabels',{''},...
                     'lineColor',[],'lineWidth',[],...
                     'xlabel','', 'ylabel','','title','');
                   
  flag_plotTitinCurves=0;

  %%
  % Data
  %%
  fileTrombitas1998Figure5 = 'experiments/TrombitasGreaserFrenchGranzier1998/Trombitas1998_Figure5.csv';
  dataTrombitas1998Figure5 = loadDigitizedData(fileTrombitas1998Figure5,...
                        'Sarcomere Length','PEVK Width (um)',...
                        {'a','b'},'Trombitas 1998 Figure 5');  
                      
  fileTrombitas1998Figure6 = 'experiments/TrombitasGreaserFrenchGranzier1998/Trombitas1998_Figure6.csv';
  dataTrombitas1998Figure6 = loadDigitizedData(fileTrombitas1998Figure6,...
                        'Sarcomere Length','PEVK Width (um)',...
                        {'a'},'Trombitas 1998 Figure 6');  

  %%
  % Tendon force-length
  %%
                   
  idxTendon =1;
  idx       =idxTendon;
  
  plotProps(idx).xlim = ...
    [1.0, (1+felineSoleusMusculotendonPropertiesUpd.tendonStrainAtOneNormForce)] ...
    +[-0.001,0.001];
  
  plotProps(idx).ylim = [0,1.0]+[-0.01,0.01];
  plotProps(idx).domain      = plotProps(idx).xlim + [-0.01,0.01];
  plotProps(idx).xticks      = [round(1.0,3), round((1+felineSoleusMusculotendonPropertiesUpd.tendonStrainAtOneNormForce),3)];
  plotProps(idx).xticklabels = {'0','$$e^T_o$$'};
  plotProps(idx).yticks      = [round(0.0,3), round(2/3,3),round(1.0,3)];
  plotProps(idx).yticklabels = {'0','$$f^T_{toe}$$','$$f^M_o$$'};
  plotProps(idx).lineColor   = [0,0,0];
  plotProps(idx).lineWidth   = 0.5;
  plotProps(idx).xlabel = 'Strain ($$\ell/\ell^{T}_{S}-1$$)';
  plotProps(idx).ylabel = 'Normalized Force ($$f/f^{M}_{o}$$)';
  plotProps(idx).title  = 'A. Tendon-force length ($$\mathbf{f}^T$$)';

  %%
  % CE force-length
  %%
  idxCELength = 2;
  idx         =idxCELength;
  normMyosinLength = felineSoleusSarcomerePropertiesUpd.normMyosinHalfLength*2;
  normMyosinBareLength = felineSoleusSarcomerePropertiesUpd.normMyosinBareHalfLength*2;
  normActinLength = felineSoleusSarcomerePropertiesUpd.normActinLength;
  normZLineThickness = felineSoleusSarcomerePropertiesUpd.normZLineLength;
  normSarcomereLengthZeroForce = felineSoleusSarcomerePropertiesUpd.normSarcomereLengthZeroForce;  
  
  %from createFiberActiveForceLength.m
  plotProps(idx).xlim = ...
    felineSoleusNormMuscleCurves.activeForceLengthCurve.xEnd ...
    +[-0.02,0.02];
  
  plotProps(idx).ylim = [0,1.0]+[-0.01,0.01];
  plotProps(idx).domain      = plotProps(idx).xlim + [-0.01,0.01];
  plotProps(idx).xticks      = [round(activeForceLengthCurveAnnotationPoints.x(1,1),2),...                                
                                1,...
                                round(felineSoleusNormMuscleCurves.fiberForceLengthCurve.xEnd(1,2),2),...
                                round(activeForceLengthCurveAnnotationPoints.x(end,1),2)];
  plotProps(idx).xticklabels = {num2str(plotProps(idx).xticks(1,1)),...
                                num2str(plotProps(idx).xticks(1,2)),...
                                num2str(plotProps(idx).xticks(1,3)),...
                                num2str(plotProps(idx).xticks(1,4))};
                              
  plotProps(idx).yticks      = [0,1];
  plotProps(idx).yticklabels = {'0','$$f^M_o$$'};
  plotProps(idx).lineColor   = [0,0,0];
  plotProps(idx).lineWidth   = 0.5;
  plotProps(idx).xlabel = 'Normalized Length ($$\ell/\ell^{M}_{o}$$)';
  plotProps(idx).ylabel = 'Normalized Force ($$f/f^{M}_{o}$$)';
  plotProps(idx).title  = 'B. CE force-length ($$\mathbf{f}^{L}\,\&\,\mathbf{f}^{PE}$$)';

  %%
  % CE force-velocity
  idxCEVelocity = 3;
  idx         =idxCEVelocity;

  %from createFiberActiveForceLength.m
  plotProps(idx).xlim = [-1.01,1.01];
  
  plotProps(idx).ylim = [0,felineSoleusNormMuscleCurves.fiberForceVelocityCurve.yEnd(1,2)] ...
                        +[-0.01,0.01];
                      
  plotProps(idx).domain      = plotProps(idx).xlim + [-0.01,0.01];
  
  plotProps(idx).xticks      = [-1,0,1];
  plotProps(idx).xticklabels = {'$$-1$$','0','$$1$$'};
                              
  plotProps(idx).yticks      = [0,1,felineSoleusNormMuscleCurves.fiberForceVelocityCurve.yEnd(1,2)];
  plotProps(idx).yticklabels = {0,1,num2str(round(plotProps(idx).yticks(1,end),2))};
  plotProps(idx).lineColor   = [0,0,0];
  plotProps(idx).lineWidth   = 0.5;
  plotProps(idx).xlabel = 'Normalized Velocity ($$v/v^{M}_{max}$$)';
  plotProps(idx).ylabel = 'Normalized Force ($$f/f^{M}_{o}$$)';
  plotProps(idx).title  = 'C. CE force-velocity ($$\mathbf{f}^{V}$$)';
  %%

  
  %%
  % PE: ECM + titin force-length
  %%
  idxPeEcmTitinLength = 4;
  idx         =idxPeEcmTitinLength;
  

  %from createFiberActiveForceLength.m
  plotProps(idx).xlim = ...
    felineSoleusNormMuscleCurves.fiberForceLengthCurve.xEnd ...
    +[-0.02,0.02];
  
  plotProps(idx).ylim = [0,1.0]+[-0.01,0.01];
  plotProps(idx).domain      = plotProps(idx).xlim + [-0.01,0.01];
  plotProps(idx).xticks      = [round(felineSoleusNormMuscleCurves.fiberForceLengthCurve.xEnd(1,1),2),...                                
                                1,...
                                round(felineSoleusNormMuscleCurves.fiberForceLengthCurve.xEnd(1,2),2)];
  plotProps(idx).xticklabels = {num2str(plotProps(idx).xticks(1,1)),...
                                num2str(plotProps(idx).xticks(1,2)),...
                                num2str(plotProps(idx).xticks(1,3))};
                              
  plotProps(idx).yticks      = [0,1];
  plotProps(idx).yticklabels = {'0','$$f^M_o$$'};
  plotProps(idx).lineColor   = [0,0,0];
  plotProps(idx).lineWidth   = 0.5;
  plotProps(idx).xlabel = 'Normalized Length ($$\ell/\ell^{M}_{o}$$)';
  plotProps(idx).ylabel = 'Normalized Force ($$f/f^{M}_{o}$$)';
  plotProps(idx).title  = 'A. ECM ($$\mathbf{f}^{ECM}$$) \& Titin ($$\mathbf{f}^{IgP}+\mathbf{f}^{PEVK}$$)';  
  
  %%
  % PE: IgP, PEVK, Titin
  %%
  idxIgpPevkForceLength = 5;
  idx         =idxIgpPevkForceLength;
  

  %from createFiberActiveForceLength.m

  
  lambdaECM = felineSoleusSarcomerePropertiesUpd.extraCellularMatrixPassiveForceFraction;  
  
  
  ligpOpt  = felineSoleusSarcomerePropertiesUpd.IGPNormLengthAtOptimalFiberLength;
  lpevkOpt = felineSoleusSarcomerePropertiesUpd.PEVKNormLengthAtOptimalFiberLength;
  
  ligpZero  = felineSoleusNormMuscleCurves.forceLengthIgpCurve.xEnd(1,1);
  lpevkZero = felineSoleusNormMuscleCurves.forceLengthPevkIgdCurve.xEnd(1,1);
  
  ligpOpt  = felineSoleusSarcomerePropertiesUpd.IGPNormLengthAtOptimalFiberLength;
  lpevkOpt = felineSoleusSarcomerePropertiesUpd.PEVKNormLengthAtOptimalFiberLength;  
  
  ligpFiso  = calcBezierYFcnXDerivative((1-lambdaECM),...
              felineSoleusNormMuscleCurves.('forceLengthIgpInverseCurve'),...
              0);
            
  lpevkFiso = calcBezierYFcnXDerivative((1-lambdaECM),...
              felineSoleusNormMuscleCurves.('forceLengthPevkIgdInverseCurve'),...
              0);

  
  normLengthZToT12 = ...
    felineSoleusSarcomerePropertiesUpd.ZLineToT12NormLengthAtOptimalFiberLength;   
  
  normLengthIgdFixed = ...
    felineSoleusSarcomerePropertiesUpd.IGDFixedNormLengthAtOptimalFiberLength;    
  
  ltitinOpt = normLengthZToT12+normLengthIgdFixed + ligpOpt+lpevkOpt;
  ltitinFiso = normLengthZToT12+normLengthIgdFixed + ligpFiso+lpevkFiso;
    
  xposXTickIgpPevkFL = [ligpOpt,ligpFiso,lpevkFiso,ltitinOpt].*2;
  textXTickIgpPevkFL = {'$$\ell^{Igp}_{o}$$','$$\ell^{PEVK}_{o}$$',...
                   '$$\ell^{Igp}_{1-\lambda}$$','$$\ell^{PEVK}_{1-\lambda}$$',...
                   '$$\ell^{T12}+\ell^{Igp}(f)+\ell^{PEVK}(f)+\ell^{Igd}$$'};  

  plotProps(idx).xlim = [0,max(ligpFiso,lpevkFiso)] + [-0.01,0.01];  
  plotProps(idx).ylim = [0,(1-lambdaECM)] + [-0.01,0.01]*(1-lambdaECM);
  plotProps(idx).domain      = plotProps(idx).xlim + [-0.01,0.01];
                 
                 
  plotProps(idx).xticks      = [0,ligpFiso, lpevkFiso];
  plotProps(idx).xticklabels = {num2str(round(plotProps(idx).xticks(1,1),2)),...
                                num2str(round(plotProps(idx).xticks(1,2),2)),...
                                num2str(round(plotProps(idx).xticks(1,3),2))};
                              
  
                              
  plotProps(idx).yticks      = [0,(1-lambdaECM)];
  plotProps(idx).yticklabels = {'0',['$$',num2str(round(1-lambdaECM,2)),'f^M_o$$']};
  plotProps(idx).lineColor   = [0,0,0];
  plotProps(idx).lineWidth   = 0.5;
  plotProps(idx).xlabel = 'Normalized Length ($$\ell/\ell^{M}_{o}$$)';
  plotProps(idx).ylabel = 'Normalized Force ($$f/f^{M}_{o}$$)';
  plotProps(idx).title  = 'B. Titin segments ($$\mathbf{f}^{IgP}$$ and $$\mathbf{f}^{PEVK}$$)';  
    
  

  %%
  % PE: ECM + titin force-length
  %%
  idxTrombitasFigure5 = 6;
  idx   =idxTrombitasFigure5;
  
  loptHuman=2.725;%um
  %from createFiberActiveForceLength.m
  plotProps(idx).xlim = [2.25,4.75];  
  plotProps(idx).ylim = [0,1.5];
  
  plotProps(idx).domain      = plotProps(idx).xlim + [-0.01,0.01];
  plotProps(idx).xticks      = [2.25:0.5:4.75];
  plotProps(idx).xticklabels = cell(1,length(plotProps(idx).xticks),1);
  for z=1:1:length(plotProps(idx).xticks)
    plotProps(idx).xticklabels{z} = num2str(round(plotProps(idx).xticks(1,z),3));
  end
  plotProps(idx).yticks      = [0:0.2:1.2];
  plotProps(idx).yticklabels = cell(1,length(plotProps(idx).yticks));
  for z=1:1:length(plotProps(idx).yticks)
    plotProps(idx).yticklabels{z} = num2str(round(plotProps(idx).yticks(1,z),3));
  end
  
  plotProps(idx).lineColor   = [0,0,0];
  plotProps(idx).lineWidth   = 0.5;
  plotProps(idx).xlabel = 'Sarcomere Length ($$\mu m$$)';
  plotProps(idx).ylabel = 'epitope to Z-line ($$\mu m$$)';
  plotProps(idx).title  = {'C. Human soleus titin segment elongation','(Trombitas et al. 1998 Fig. 5)'};  
    
  
  
  
  %%
  %Tendon
  %%
  figure(fig_pubCurves);
  subplotTendon = reshape(subPlotPanel(1,1,:),1,4);
  subplot('Position',subplotTendon);
  
    idx=idxTendon;    
    curveSample = calcBezierYFcnXCurveSampleVector(...
                  felineSoleusNormMuscleCurves.('tendonForceLengthCurve'), ...
                  200,plotProps(idx).domain);

    plot(curveSample.x, ...
      curveSample.y,...
      '-','Color',plotProps(idx).lineColor,...
          'LineWidth',plotProps(idx).lineWidth);
    hold on;

    xticks(plotProps(idx).xticks);
    xticklabels(plotProps(idx).xticklabels);

    yticks(plotProps(idx).yticks);
    yticklabels(plotProps(idx).yticklabels);

    text( felineSoleusNormMuscleCurves.tendonForceLengthCurve.xEnd(1,2),...
          felineSoleusNormMuscleCurves.tendonForceLengthCurve.yEnd(1,2),...
          '$$k^{T}_o$$',...
          'VerticalAlignment','top',...
          'HorizontalAlignment','left',...
          'FontSize',8);

    hold on;

    plot(felineSoleusNormMuscleCurves.tendonForceLengthCurve.xEnd(1,1),...
          felineSoleusNormMuscleCurves.tendonForceLengthCurve.yEnd(1,1),...
          '.','Color',plotProps(idx).lineColor);
    hold on;    
    
    plot(felineSoleusNormMuscleCurves.tendonForceLengthCurve.xEnd(1,2),...
          felineSoleusNormMuscleCurves.tendonForceLengthCurve.yEnd(1,2),...
          '.','Color',plotProps(idx).lineColor);
    hold on;

    xlim(plotProps(idx).xlim);
    ylim(plotProps(idx).ylim);

    box off;

    xlabel(plotProps(idx).xlabel);
    ylabel(plotProps(idx).ylabel);
    title(plotProps(idx).title);
    
  %%
  %CE force-length
  %%
  figure(fig_pubCurves);  
  subplotCELength = reshape(subPlotPanel(1,2,:),1,4);
  subplot('Position',subplotCELength);  
  
    idx=idxCELength;    
    curveSampleFL = calcBezierYFcnXCurveSampleVector(...
                  felineSoleusNormMuscleCurves.('activeForceLengthCurve'), ...
                  200,plotProps(idx).domain);

    plot(curveSampleFL.x, ...
      curveSampleFL.y,...
      '-','Color',plotProps(idx).lineColor,...
          'LineWidth',plotProps(idx).lineWidth);
    hold on;

    curveSampleFPE = calcBezierYFcnXCurveSampleVector(...
                  felineSoleusNormMuscleCurves.('fiberForceLengthCurve'), ...
                  200,plotProps(idx).domain);

    plot(curveSampleFPE.x, ...
      curveSampleFPE.y,...
      '-','Color',[1,1,1],...
          'LineWidth',plotProps(idx).lineWidth*2);
    hold on;    
    
    plot(curveSampleFPE.x, ...
      curveSampleFPE.y,...
      '--','Color',plotProps(idx).lineColor,...
          'LineWidth',plotProps(idx).lineWidth);
    hold on;    
        
    plot(felineSoleusActiveForceLengthDataUpd(:,1),...
         felineSoleusActiveForceLengthDataUpd(:,2),...
         '+','Color',[0,0,1],'MarkerSize',5,...
         'LineWidth',plotProps(idx).lineWidth);
    hold on;

    plot(felineSoleusPassiveForceLengthDataUpd(:,1),...
         felineSoleusPassiveForceLengthDataUpd(:,2),...
         'x','Color',[0,0,1],'MarkerSize',5,...
         'LineWidth',plotProps(idx).lineWidth);
    hold on;
    
     
%     plot(felineSoleusActiveForceLengthLineBestFit(:,1),...
%          felineSoleusActiveForceLengthLineBestFit(:,2),...
%          '-','Color',[1,0,0],...
%          'LineWidth',plotProps(idx).lineWidth);
%     hold on;
    
    xticks(plotProps(idx).xticks);
    xticklabels(plotProps(idx).xticklabels);

    yticks(plotProps(idx).yticks);
    yticklabels(plotProps(idx).yticklabels);

    xlim(plotProps(idx).xlim);
    ylim(plotProps(idx).ylim);

    text( plotProps(idx).xticks(1,3),0.6,...
          'Herzog 2002',...
          'FontSize',8,...
          'Color',[0,0,1]);

    
    box off;

    xlabel(plotProps(idx).xlabel);
    ylabel(plotProps(idx).ylabel);
    title(plotProps(idx).title);    
  
  %%
  % Force-velocity
  %%
  figure(fig_pubCurves);  
  subplotCEVelocity = reshape(subPlotPanel(1,3,:),1,4);
  subplot('Position',subplotCEVelocity);  
  
    idx=idxCEVelocity;    
    curveSampleFV = calcBezierYFcnXCurveSampleVector(...
                  felineSoleusNormMuscleCurves.('fiberForceVelocityCurve'), ...
                  200,plotProps(idx).domain);

    plot(curveSampleFV.x, ...
      curveSampleFV.y,...
      '-','Color',plotProps(idx).lineColor,...
          'LineWidth',plotProps(idx).lineWidth);
    hold on;        
    
    plot(0,...
         1,...
         '.','Color',plotProps(idx).lineColor);
    hold on;    
    
    text(0.01,-0.01, 'Lengthening','FontSize',8,...
         'VerticalAlignment','bottom');
    hold on;
    text(-0.01,-0.01, 'Shortening','FontSize',8,...
         'HorizontalAlignment','Right',...
         'VerticalAlignment','bottom');
    hold on;
    
    xticks(plotProps(idx).xticks);
    xticklabels(plotProps(idx).xticklabels);

    yticks(plotProps(idx).yticks);
    yticklabels(plotProps(idx).yticklabels);

    xlim(plotProps(idx).xlim);
    ylim(plotProps(idx).ylim);

    box off;

    xlabel(plotProps(idx).xlabel);
    ylabel(plotProps(idx).ylabel);
    title(plotProps(idx).title);    
    
if(flag_plotTitinCurves==1)    
  %%
  % Passive curve decomposition: total, ecm, titin
  %%    
  figure(fig_pubCurves);  
  subplotPeEcmTitin = reshape(subPlotPanel(2,1,:),1,4);
  subplot('Position',subplotPeEcmTitin);  
  idx=idxPeEcmTitinLength;
  
  fpeDomain = felineSoleusNormMuscleCurves.('fiberForceLengthCurve').xEnd ...
             +[-0.01,0.01]+[0,1];
  curveSampleFPE = calcBezierYFcnXCurveSampleVector(...
                felineSoleusNormMuscleCurves.('fiberForceLengthCurve'), ...
                200,fpeDomain);

  %Now go and evaluate the ECM and the titin elements at all of the equivalent
  %lengths
  z0 = zeros(size(curveSampleFPE.x));
  curveSampleECMHalf  = struct('x',z0,'y',z0,'dydx',z0,'d2ydx2',z0,...
                           'd3ydx3',z0,'intYdx',[]);
  curveSampleIgp  = struct('x',z0,'y',z0,'dydx',z0,'d2ydx2',z0,...
                           'd3ydx3',z0,'intYdx',[]);
  curveSamplePevkIgd = struct('x',z0,'y',z0,'dydx',z0,'d2ydx2',z0,...
                           'd3ydx3',z0,'intYdx',[]);
  curveSampleTitin = struct('x',z0,'y',z0,'dydx',z0,'d2ydx2',z0,...
                           'd3ydx3',z0,'intYdx',[]);
  curveSampleTitinActive = struct('x',z0,'y',z0,'dydx',z0,'d2ydx2',z0,...
                           'd3ydx3',z0,'intYdx',[]);  
  curveSampleIgpB  = struct('x',z0,'y',z0,'dydx',z0,'d2ydx2',z0,...
                           'd3ydx3',z0,'intYdx',[]);
  curveSamplePevkB = struct('x',z0,'y',z0,'dydx',z0,'d2ydx2',z0,...
                           'd3ydx3',z0,'intYdx',[]);
  curveSampleTitinB = struct('x',z0,'y',z0,'dydx',z0,'d2ydx2',z0,...
                           'd3ydx3',z0,'intYdx',[]);
  curveSampleTitinActiveB = struct('x',z0,'y',z0,'dydx',z0,'d2ydx2',z0,...
                           'd3ydx3',z0,'intYdx',[]);                           
  normLengthZToT12 = ...
    felineSoleusSarcomerePropertiesUpd.ZLineToT12NormLengthAtOptimalFiberLength;   
  
  normLengthIgpAtOptimalFiberLength = ...
    felineSoleusSarcomerePropertiesUpd.IGPNormLengthAtOptimalFiberLength;
  
  normLengthIgdFixedHuman = ...
    felineSoleusSarcomerePropertiesUpd.IGDFixedNormLengthAtOptimalFiberLengthHuman; 

  normLengthIgdFixed = ...
    felineSoleusSarcomerePropertiesUpd.IGDFixedNormLengthAtOptimalFiberLength; 
  
  
  normLengthIgdTotal = ...
    felineSoleusSarcomerePropertiesUpd.IGDTotalNormLengthAtOptimalFiberLength;
  
  lengthZ2Igp  = z0;
  lengthZ2Pevk = z0;
  lengthZ2Igd  = z0;                  
  lengthCE     = z0;
  
  lengthZ2IgpB  = z0;
  lengthZ2PevkB = z0;
  lengthZ2IgdB  = z0;                  
  lengthCEB     = z0;
  
  
[ curveSampleECMHalf,...
  curveSampleTitin,...
  curveSampleTitinActive,...
  curveSampleIgp,...
  curveSamplePevkIgd] = ...
  sampleTitinCurves(curveSampleFPE.x.*0.5,...
                    felineSoleusNormMuscleCurves.forceLengthECMHalfCurve,...
                    felineSoleusNormMuscleCurves.forceLengthIgpCurve,...
                    felineSoleusNormMuscleCurves.forceLengthIgpInverseCurve,...
                    felineSoleusNormMuscleCurves.forceLengthPevkIgdCurve,...
                    felineSoleusNormMuscleCurves.forceLengthPevkIgdInverseCurve,...
                    normLengthZToT12,...
                    normLengthIgpAtOptimalFiberLength,...
                    normLengthIgdFixed); 

  lengthZ2Igp  = normLengthZToT12.*ones(size(curveSampleFPE.x));
  lengthZ2Pevk = lengthZ2Igp + curveSampleIgp.x;
  lengthZ2Igd  = lengthZ2Pevk + curveSamplePevkIgd.x;                  
  lengthCE     = curveSampleFPE.x.*0.5;                  
                  
[ curveSampleECMHalfHuman,...
  curveSampleTitinHuman,...
  curveSampleTitinActiveHuman,...
  curveSampleIgpHuman,...
  curveSamplePevkIgdHuman] = ...
  sampleTitinCurves(curveSampleFPE.x.*0.5,...
                    felineSoleusNormMuscleCurves.forceLengthECMHalfCurve,...
                    felineSoleusNormMuscleCurves.forceLengthIgpCurveHuman,...
                    felineSoleusNormMuscleCurves.forceLengthIgpInverseCurveHuman,...
                    felineSoleusNormMuscleCurves.forceLengthPevkCurveHuman,...
                    felineSoleusNormMuscleCurves.forceLengthPevkInverseCurveHuman,...
                    normLengthZToT12,...
                    normLengthIgpAtOptimalFiberLength,...
                    normLengthIgdFixedHuman);    
                  
  lengthZ2IgpHuman  = normLengthZToT12.*ones(size(curveSampleFPE.x));
  lengthZ2PevkHuman = lengthZ2IgpHuman + curveSampleIgpHuman.x;
  lengthZ2IgdHuman  = lengthZ2PevkHuman + curveSamplePevkIgdHuman.x;                  
  lengthCEHuman     = curveSampleFPE.x.*0.5;                  
                  
[ curveSampleECMHalfHumanB,...
  curveSampleTitinHumanB,...
  curveSampleTitinActiveHumanB,...
  curveSampleIgpHumanB,...
  curveSamplePevkIgdHumanB] = ...
  sampleTitinCurves(curveSampleFPE.x.*0.5,...
                    felineSoleusNormMuscleCurves.forceLengthECMHalfCurve,...
                    felineSoleusNormMuscleCurves.forceLengthIgpCurveHumanB,...
                    felineSoleusNormMuscleCurves.forceLengthIgpInverseCurveHumanB,...
                    felineSoleusNormMuscleCurves.forceLengthPevkCurveHumanB,...
                    felineSoleusNormMuscleCurves.forceLengthPevkInverseCurveHumanB,...
                    normLengthZToT12,...
                    normLengthIgpAtOptimalFiberLength,...
                    normLengthIgdTotal);    
                  
  lengthZ2IgpHumanB  = normLengthZToT12.*ones(size(curveSampleFPE.x));
  lengthZ2PevkHumanB = lengthZ2IgpHumanB + curveSampleIgpHumanB.x;
  lengthZ2IgdHumanB  = lengthZ2PevkHumanB + curveSamplePevkIgdHumanB.x;                  
  lengthCEHumanB     = curveSampleFPE.x.*0.5;       
  
              
              

  fill([curveSampleECMHalf.x;fliplr(curveSampleECMHalf.x')'].*2, ...
       [zeros(size(curveSampleECMHalf.y)); fliplr(curveSampleECMHalf.y')'],...
        [0.9,0.9,1],'EdgeColor','none');
  hold on;  

  plot( curveSampleECMHalf.x.*2, ...
        curveSampleECMHalf.y,...
        '-','Color',[0,0,1],...
        'LineWidth',plotProps(idx).lineWidth*2);
  hold on;    
 
  
  colorTitin = [1,0,1];
  colorTitinLight = colorTitin.*0.1 + [1,1,1].*0.9;
  
  fill([curveSampleTitin.x;fliplr(curveSampleTitin.x')'].*2, ...
       [curveSampleECMHalf.y; fliplr(curveSampleTitin.y'+curveSampleECMHalf.y')'],...
        colorTitinLight,'EdgeColor','none');
  hold on;  

  plot( curveSampleTitin.x.*2, ...
        curveSampleTitin.y+curveSampleECMHalf.y,...
        '-','Color',colorTitin,...
        'LineWidth',plotProps(idx).lineWidth);
  hold on;     
  
  
  plot( curveSampleFPE.x, ...
        curveSampleFPE.y,...
        '--','Color',[0,0,0],...
        'LineWidth',plotProps(idx).lineWidth);
  hold on;    
  
  plot( curveSampleTitinActive.x.*2,...
        curveSampleTitinActive.y+curveSampleECMHalf.y,...
        '-','Color',colorTitin,...
        'LineWidth',plotProps(idx).lineWidth);
  hold on;   
  
  yTxt = 1;
  xTxt = interp1(curveSampleTitinActive.y+curveSampleECMHalf.y,...
                 curveSampleTitinActive.x*2, yTxt);
  hTa =  text(xTxt,yTxt,'Titin (active)',...
       'FontSize',8,...
       'HorizontalAlignment','right',...
       'VerticalAlignment','bottom');
  hold on
  angle = atan(curveSampleTitinActive.dydx(end)...
               +curveSampleECMHalf.dydx(end))*(180/pi);  
             
  set(hTa,'Rotation',angle+labelRotationOffset);
  
  lambdaECM = felineSoleusSarcomerePropertiesUpd.extraCellularMatrixPassiveForceFraction;
  
  text(1.36,...
       0,...
       ['$$\lambda^{ECM}=',sprintf('%1.2f',lambdaECM),'$$'],...
       'FontSize',8,...
       'HorizontalAlignment','right',...
       'VerticalAlignment','bottom');
  hold on;
  
  plot(felineSoleusPassiveForceLengthDataUpd(:,1),...
       felineSoleusPassiveForceLengthDataUpd(:,2),...
       'x','Color',[0,0,1],'MarkerSize',5,...
       'LineWidth',plotProps(idx).lineWidth);
  hold on;  
  
  
  yTxt = lambdaECM;
  
  yErr = Inf;
  z=length(curveSampleECMHalf.y);
  while z > 2 && ...
      ((curveSampleECMHalf.y(z,1)-yTxt)*(curveSampleECMHalf.y(z-1,1)-yTxt) > 0)
    z=z-1;
  end
  
  xTxt = curveSampleECMHalf.x(z,1);
  yTxt = curveSampleECMHalf.y(z,1);
  angle = atan(curveSampleECMHalf.dydx(z,1))*(180/pi);
  
  hecm= text(xTxt*2, yTxt,...
       'ECM          ',...
       'FontSize',8,...
       'HorizontalAlignment','right',...
       'VerticalAlignment','top');
  hold on;
  set(hecm,'Rotation',angle+labelRotationOffset*2.5);
  hold on;
  
  hti = text(xTxt*2, yTxt,...
     'Titin (passive)',...
     'FontSize',8,...
     'HorizontalAlignment','right',...
     'VerticalAlignment','bottom');
  hold on;
  
  set(hti,'Rotation',angle+labelRotationOffset*2.5);
  hold on;
    
  box off

  
  
  
  xticks(plotProps(idx).xticks);
  xticklabels(plotProps(idx).xticklabels);

  yticks(plotProps(idx).yticks);
  yticklabels(plotProps(idx).yticklabels);

  xlim(plotProps(idx).xlim);
  ylim(plotProps(idx).ylim);

  box off;

  xlabel(plotProps(idx).xlabel);
  ylabel(plotProps(idx).ylabel);
  title(plotProps(idx).title);      
  
  %%
  % Igp + PEVK force-length
  %%    
  figure(fig_pubCurves);  
  subplotIgpPevk = reshape(subPlotPanel(2,2,:),1,4);
  subplot('Position',subplotIgpPevk);  
  idx=idxIgpPevkForceLength;
    
  %textXTickIgpPevkFL
  
  colorIgp   = [1,0,0];
  colorPEVK  = [0.5,0,0];
  
  plot( curveSampleIgp.x,...
        curveSampleIgp.y,...
        '-','Color',colorIgp,...
        'LineWidth',plotProps(idx).lineWidth);
  hold on;
  
  angleStretch  = 1.1;
  
  eIsoIgp = ligpFiso/ligpOpt;

  yTxt = (1-lambdaECM)*0.9;
  xTxt = interp1( curveSampleIgp.y, curveSampleIgp.x, yTxt);
  
  hIgp = text(xTxt, yTxt,...
        ['$$\mathbf{f}^{IgP}$$'],'FontSize',8,...
        'HorizontalAlignment','right',...
        'VerticalAlignment','bottom');     
  hold on;

  yTxt = (1-lambdaECM)*0.7;
  xTxt = interp1( curveSampleIgp.y, curveSampleIgp.x, yTxt);  
  
  kigp = felineSoleusNormMuscleCurves.forceLengthIgpCurve.dydxEnd(1,2);
  higp = text(xTxt,yTxt,...
       ['$$\hat{k}^{IgP}_o=',num2str(round(kigp,2)),'$$'],'FontSize',8,...
       'HorizontalAlignment','right',...
       'VerticalAlignment','bottom');      
  hold on;   
  
  %angleIgp = atan(curveSampleIgp.dydx(end))*(180/pi);
  %set(hIgp,'Rotation',angleIgp+labelRotationOffset);
  %hIgp2 = text(curveSampleIgp.x(end).*2,...
  %      curveSampleIgp.y(end),...
  %      ['$$e_{ISO}=\, ',num2str(round(eIsoIgp,2)),'$$'],'FontSize',8,...
  %      'HorizontalAlignment','right',...
  %      'VerticalAlignment','top');           
  %angleIgp2 = atan(curveSampleIgp.dydx(end))*(180/pi);
  %set(hIgp2,'Rotation',angleIgp2*angleStretch);
  
  
  plot( curveSamplePevkIgd.x,...
        curveSamplePevkIgd.y,...
        '-','Color',colorPEVK,...
        'LineWidth',plotProps(idx).lineWidth);
  hold on;  

  yTxt = (1-lambdaECM)*0.9;
  xTxt = interp1( curveSamplePevkIgd.y, curveSamplePevkIgd.x, yTxt);  
  
  %eIsoPevk = lpevkFiso/lpevkOpt;
  hpevk = text(xTxt,yTxt,...
       ['$$\mathbf{f}^{PEVK}$$'],'FontSize',8,...
       'HorizontalAlignment','right',...
       'VerticalAlignment','bottom');      
  hold on;  
  
  yTxt = (1-lambdaECM)*0.7;
  xTxt = interp1( curveSamplePevkIgd.y, curveSamplePevkIgd.x, yTxt);  
  
  kpevk = felineSoleusNormMuscleCurves.forceLengthPevkIgdCurve.dydxEnd(1,2);
  hpevk = text(xTxt,yTxt,...
       ['$$\hat{k}^{PEVK}_o=',num2str(round(kpevk,2)),'$$'],'FontSize',8,...
       'HorizontalAlignment','right',...
       'VerticalAlignment','bottom');      
  hold on;    
  
  %anglePevk = atan(curveSamplePevkIgd.dydx(end))*(180/pi);
  %set(hpevk,'Rotation',anglePevk+labelRotationOffset);  

  %hpevk2 = text(curveSamplePevkIgd.x(end).*2,...
  %      curveSamplePevkIgd.y(end),...
  %      ['$$e_{ISO}=\,',num2str(round(eIsoPevk,2)),' $$'],'FontSize',8,...
  %      'HorizontalAlignment','right',...
  %      'VerticalAlignment','top');      
  %anglePevk = atan(curveSamplePevkIgd.dydx(end))*(180/pi);
  %set(hpevk2,'Rotation',anglePevk*angleStretch);  
  %hold on;
  
%   plot( curveSampleTitin.x.*2,...
%         curveSampleTitin.y,...
%         '-','Color',colorTitin,...
%         'LineWidth',plotProps(idx).lineWidth);
%   hold on;
%   
%   eIsoTitin = (ligpFiso+lpevkFiso+normLengthZToT12 + normLengthIgdFixed)...
%              /(ligpOpt+lpevkOpt+normLengthZToT12 + normLengthIgdFixed);
% 
%   angleStretchBig=1.2;
% 
%   yTxt = 1-lambdaECM;
%   xTxt = interp1( curveSampleTitin.y, curveSampleTitin.x, yTxt);    
%   
%   hti = text(xTxt.*2,...
%         yTxt,...
%         'Titin (passive)',...
%         'FontSize',8,...
%         'HorizontalAlignment','right',...
%         'VerticalAlignment','bottom');
%   angleTi = atan(curveSampleTitin.dydx(end))*(180/pi);
%   set(hti,'Rotation',angleTi*angleStretchBig);
%   
%   hold on;
  
%   plot( curveSampleTitinActive.x.*2,...
%         curveSampleTitinActive.y,...
%         '--','Color',colorTitin,...
%         'LineWidth',plotProps(idx).lineWidth);
%   hold on;  
% 
% 
%   yTxt = 1-lambdaECM;
%   xTxt = interp1( curveSampleTitinActive.y, curveSampleTitinActive.x, yTxt); 
%   
%   hti = text(xTxt.*2,...
%         yTxt,...
%         'Titin (active)',...
%         'FontSize',8,...
%         'HorizontalAlignment','right',...
%         'VerticalAlignment','bottom');
%   angleTi = atan(curveSampleTitinActive.dydx(end))*(180/pi);
%   set(hti,'Rotation',angleTi*angleStretch);
  
%   angleStretchBig=1.2;
%   hti = text(curveSampleTitin.x(end).*2,...
%         curveSampleTitin.y(end),...
%         '$$\ell^{T12}+\ell^{Igp}(f)+\ell^{PEVK}(f)+\ell^{Igd}$$',...
%         'FontSize',8,...
%         'HorizontalAlignment','right',...
%         'VerticalAlignment','bottom');
%   angleTi = atan(curveSampleTitin.dydx(end))*(180/pi);
%   set(hti,'Rotation',angleTi*angleStretchBig);
  
%   hti2=text(curveSampleTitin.x(end).*2,...
%         curveSampleTitin.y(end),...
%         ['$$e_{ISO}=',num2str(round(eIsoTitin,2)),'$$'],...
%         'FontSize',8,...
%         'HorizontalAlignment','right',...
%         'VerticalAlignment','top');
%   angleTi2 = atan(curveSampleTitin.dydx(end))*(180/pi);
%   set(hti2,'Rotation',angleTi2*angleStretchBig);
%   
%   
%   hold on;
  
%   text(xposXTickIgpPevkFL(1,1),...
%        0.05,...
%        textXTickIgpPevkFL{1,1},...
%        'FontSize',8')
%   hold on;
%   
%   text(xposXTickIgpPevkFL(2,1),...
%        0.05,...
%        textXTickIgpPevkFL{1,1},...
%        'FontSize',8')
%   hold on;

  
  
  %xposXTickIgpPevkFL = [ligpOpt,lpevkOpt,ligpFiso,lpevkFiso];
  %textXTickIgpPevkFL = {'$$\ell^{Igp}_{o}$$','$$\ell^{PEVK}_{o}$$',...
  %                 '$$\ell^{Igp}_{1-\lambda}$$','$$\ell^{PEVK}_{1-\lambda}$$'};  
  
  box off

  xticks(plotProps(idx).xticks);
  xticklabels(plotProps(idx).xticklabels);

  yticks(plotProps(idx).yticks);
  yticklabels(plotProps(idx).yticklabels);

  xlim(plotProps(idx).xlim);
  ylim(plotProps(idx).ylim);

  xlabel(plotProps(idx).xlabel);
  ylabel(plotProps(idx).ylabel);
  title(plotProps(idx).title);   
  
  %%
  % Igp + PEVK force-length
  %%    
  figure(fig_pubCurves);  
  subplotTrombitasFigure5 = reshape(subPlotPanel(2,3,:),1,4);
  subplot('Position',subplotTrombitasFigure5);  
  idx=idxTrombitasFigure5;  
  
  plot( dataTrombitas1998Figure5(2).x,...
        dataTrombitas1998Figure5(2).y,...
        'o','MarkerSize',3, ...
        'MarkerFaceColor',[1,1,1].*0.9,...
        'Color',[1,1,1].*0.75,...
        'DisplayName','Exp: Z-line (ZL) to N-end 9D10');
  hold on;

  plot( dataTrombitas1998Figure5(1).x,...
        dataTrombitas1998Figure5(1).y,...
        'o','MarkerSize',3, ...
        'MarkerFaceColor',[1,1,1].*0.5,...
        'Color',[1,1,1].*0.5,...
        'DisplayName','Exp: ZL to C-end 9D10');
  hold on;
  
  
  plot( (lengthCEHuman.*2).*loptHuman,...
        (lengthZ2PevkHuman).*loptHuman,...
        '-','Color',colorIgp,...
        'LineWidth',plotProps(idx).lineWidth,...
        'DisplayName','Model: ZL to IgD/PEVK');
      
  hold on;
  

        
  plot( (lengthCEHuman.*2).*loptHuman,...
        (lengthZ2IgdHuman).*loptHuman,...
        '-','Color',colorPEVK,...
        'LineWidth',plotProps(idx).lineWidth,...
        'DisplayName','Model: ZL to IgD/Myosin');
      
        
  plot( (lengthCEHuman.*2).*loptHuman,...
        (lengthZ2IgdHumanB).*loptHuman,...
        '--','Color',colorPEVK,...
        'LineWidth',plotProps(idx).lineWidth,...
        'DisplayName','Model: ZL to PEVK/IgD');
      
  hold on;

  legend('Location','NorthWest');
  legend boxoff;
  
  
  %text(3.25,0,'Trombitas et al. 1998',...
  %     'FontSize',8,...
  %     'HorizontalJustification','left',...
  %     'Vertical Justification','top');
  %hold on;
  
 %   lengthZ2Igp(i,1) = normLengthZToT12;
 %   lengthZ2Pevk(i,1) = lengthZ2IgpP(i,1) + ligpH;
 %   lengthZ2Igd(i,1) = lengthZ2Pevk(i,1) + lpevkH;                  
 %   lengthCE(i,1)     = xH;  
  box off

  xticks(plotProps(idx).xticks);
  xticklabels(plotProps(idx).xticklabels);

  yticks(plotProps(idx).yticks);
  yticklabels(plotProps(idx).yticklabels);

  xlim(plotProps(idx).xlim);
  ylim(plotProps(idx).ylim);

  xlabel(plotProps(idx).xlabel);
  ylabel(plotProps(idx).ylabel);
  title(plotProps(idx).title);   
   
  
end  
  here=1;

  figure(fig_pubCurves);  
  configPlotExporter;
  %print('-dpdf', [pubOutputFolder,'fig_Pub_MuscleCurves.pdf']); 