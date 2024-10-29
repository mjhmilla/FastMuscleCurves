function [figPub,pageWidth,pageHeight] ...
            = plotVexatMusculotendonCurvesPub(muscleModelStruct,...
                        normMuscleQuadraticCurves,...
                        activeForceLengthCurveAnnotationPoints,...
                        activeForceLengthData,...
                        passiveForceLengthData,...
                        flag_addMat156,...
                        flag_genericCurveConfig,...
                        figPub)

    figure(figPub);



    numberOfHorizontalPlotColumns   = 3;
    numberOfVerticalPlotRows        = 2;
    plotWidth           = 4.5;
    plotHeight          = 4.5;
    plotHorizMarginCm   = 1.5;
    plotVertMarginCm    = 2.0;

    [subPlotPanel,pageWidth,pageHeight]  = ...
        plotConfigGeneric(numberOfHorizontalPlotColumns, ...
                          numberOfVerticalPlotRows,...
                          plotWidth,...
                          plotHeight,...
                          plotHorizMarginCm,...
                          plotVertMarginCm);   

    npts=100;

    falSample = calcQuadraticBezierYFcnXCurveSampleVector(...
        normMuscleQuadraticCurves.activeForceLengthCurve, npts, []);

    falCalSample = calcQuadraticBezierYFcnXCurveSampleVector(...
        normMuscleQuadraticCurves.activeForceLengthCalibratedCurve, npts, []);

    fpeSample = calcQuadraticBezierYFcnXCurveSampleVector(...
        normMuscleQuadraticCurves.fiberForceLengthCurve, npts, []);
    fpeDomain = normMuscleQuadraticCurves.fiberForceLengthCurve.xEnd;

    lpeNZero= normMuscleQuadraticCurves.fiberForceLengthCurve.xEnd(1,1);
    lpeNOne = calcQuadraticBezierYFcnXDerivative(1.0, ...
                normMuscleQuadraticCurves.fiberForceLengthInverseCurve,...
                0);

    fvSample = calcQuadraticBezierYFcnXCurveSampleVector(...
        normMuscleQuadraticCurves.fiberForceVelocityCurve, npts, [-1.1,1.1]);

    fvCalSample = calcQuadraticBezierYFcnXCurveSampleVector(...
        normMuscleQuadraticCurves.fiberForceVelocityCalibratedCurve, npts, [-1.1,1.1]);

    ltNOne = calcQuadraticBezierYFcnXDerivative(1.0, ...
                normMuscleQuadraticCurves.tendonForceLengthInverseCurve,...
                0);

    ltNZero = 1.0;
    ftSample = calcQuadraticBezierYFcnXCurveSampleVector(...
        normMuscleQuadraticCurves.tendonForceLengthCurve, ...
        npts, [ltNZero-0.1*(ltNOne-ltNZero),ltNOne+(ltNOne-ltNZero)*0.1]);

    samples = [min(fpeDomain):((max(fpeDomain)-min(fpeDomain))/(npts-1)):max(fpeDomain)]';

    [ curveSampleECMHalf,...
      curveSampleTitin,...
      curveSampleTitinActive,...
      curveSampleIgp,...
      curveSamplePevk,...
      curveSampleIgd,...
      curveSampleProximalTitin,...
      curveSampleDistalTitin] = ...
      sampleTitinCurves(samples.*0.5,...
                        muscleModelStruct.curves.forceLengthECMHalfCurve,...
                        muscleModelStruct.curves.forceLengthProximalTitinCurve,...
                        muscleModelStruct.curves.forceLengthProximalTitinInverseCurve,...
                        muscleModelStruct.curves.forceLengthDistalTitinCurve,...
                        muscleModelStruct.curves.forceLengthDistalTitinInverseCurve,...                    
                        muscleModelStruct.curves.forceLengthIgPTitinCurve,...
                        muscleModelStruct.curves.forceLengthIgPTitinInverseCurve,...
                        muscleModelStruct.curves.forceLengthPevkTitinCurve,...
                        muscleModelStruct.curves.forceLengthPevkTitinInverseCurve,...
                        muscleModelStruct.curves.forceLengthIgDTitinCurve,...
                        muscleModelStruct.curves.forceLengthIgDTitinInverseCurve,...
                        muscleModelStruct.sarcomere.ZLineToT12NormLengthAtOptimalFiberLength,...
                        muscleModelStruct.sarcomere.IGPNormLengthAtOptimalFiberLength,...
                        muscleModelStruct.sarcomere.PEVKNormLengthAtOptimalFiberLength,...
                        muscleModelStruct.sarcomere.IGDFixedNormLengthAtOptimalFiberLength,...
                        muscleModelStruct.sarcomere.titinModelType);  

    ecmForceFraction = ...
        muscleModelStruct.sarcomere.extraCellularMatrixPassiveForceFraction;      

    fpeSample.x = curveSampleECMHalf.x.*2;
    fpeSample.y = curveSampleTitin.y.*(1-ecmForceFraction) ...
                       +curveSampleECMHalf.y.*(ecmForceFraction); 


    greyA = [0,0,0];
    greyB = [1,1,1].*0.5;
    
    blueA  = [0, 0.4470, 0.7410];
    blueB  = blueA.*0.5 + [1,1,1].*0.5;
    
    greenA = [0, 0.75, 0.75];
    greenB = greenA.*0.5 + [1,1,1].*0.5;
    
    maroonA= [0.6350, 0.0780, 0.1840];
    maroonB= maroonA.*0.5 + [1,1,1].*0.5;
 
    magentaA = [0.75, 0, 0.75];
    magentaB = magentaA.*0.5 + [1,1,1].*0.5;
    
    redA = [193, 39, 45]./255;
    redB = redA.*0.5 + [1,1,1].*0.5;
    blueB       = [0, 0.4470, 0.7410].*0.66 + [1,1,1].*0.33;
    blueC       = [0, 0.4470, 0.7410].*0.33 + [1,1,1].*0.66;
    magentaA    = [0.75, 0, 0.75];
    magentaBlue = magentaA.*0.5 + blueA.*0.5;

    maroonA     = [0.6350, 0.0780, 0.1840];    
    greenA      = [0, 0.75, 0.75];

    colorModel    = blueA;
    colorMAT156   = redA;
    colorData     = [0,0,0];



   
    magenta0 = ([170,51,119]./255);
    magenta1 = magenta0.*0.75 + [1,1,1].*0.25;
    magenta2 = magenta0.*0.75 + [0,0,0].*0.25;

    colorPassive=magenta0;

    colorEcm = magenta1;
    colorTitin=magenta2;

    colorTitinActive = [238,51,119]./255;

    lineWidthModel  = 1.5;
    lineWidthMAT156 = 1;
    lineWidthData   = 1;

    labelModel  = 'VEXAT ';
    labelMAT156 = 'MAT156 ';
    labelData   = 'Exp. ';

    if(flag_genericCurveConfig==1)
        labelModel  = '';
        labelMAT156 = '';
        flag_addMat156=0;
    end


    xBridgeDelta = ones(size(falCalSample.y)) ./ (muscleModelStruct.sarcomere.normCrossBridgeStiffness);
    xBridgeDelta = xBridgeDelta.*2;
    flag_addCrossbridgeStrain  = 0;
    flag_plotActiveTitinForces = 1;

    %%
    %Evaluate titin curves
    %%


    subplot('Position',reshape(subPlotPanel(1,1,:),1,4));

        %Active-force length plots
        if(flag_addCrossbridgeStrain==1)
            plot(falCalSample.x,falCalSample.y,...
                '-','Color',colorModel,...
                'LineWidth',lineWidthModel,...
                'DisplayName',[labelModel,'$$f^{L}_{CAL}$$']);
        else
            plot(falSample.x,falSample.y,...
                '-','Color',colorModel,...
                'LineWidth',lineWidthModel,...
                'DisplayName',[labelModel,'$$f^{L}$$']);
        end
        hold on;    
        if(flag_addMat156==1)
            plot(falSample.x,falSample.y,'--','Color',colorMAT156,...
                'LineWidth',lineWidthMAT156,...
                'DisplayName',[labelMAT156,'$$f^{L}_{156}$$']);
            hold on;
        end
        if(~isempty(activeForceLengthData))
            plot(activeForceLengthData(:,1),...
                 activeForceLengthData(:,2),...
                 'o','Color',colorData,...
                 'LineWidth',lineWidthData,...
                 'DisplayName',[labelData,'$$f^{L}$$'],...
                 'MarkerSize',4,...
                 'MarkerFaceColor',[1,1,1]);
        end


        plot(fpeSample.x,fpeSample.y,'-','Color',colorPassive,...
            'LineWidth',lineWidthModel,...
            'DisplayName',[labelMAT156,'$$f^{PE}$$']);
        hold on;


        if(~isempty(passiveForceLengthData))
            plot(passiveForceLengthData(:,1),...
                 passiveForceLengthData(:,2),...
                 's','Color',colorData,...
                 'LineWidth',lineWidthData,...
                 'DisplayName',[labelData,': $$f^{PE}$$'],...
                 'MarkerSize',4,...
                 'MarkerFaceColor',[1,1,1]);    
        end






        title('A. Active \& Passive Force-Length Relations');
        
        if(flag_genericCurveConfig==1)
            xtickSet = sort(unique(round(activeForceLengthCurveAnnotationPoints.x([1,3,4],1)',1)));
            xtickLabelSet={''};
            for k=1:1:length(xtickSet)
                tickLabel = sprintf('%1.1f%s',xtickSet(1,k),'$$\ell_o^M$$');
                if(k==1)
                    xtickLabelSet={tickLabel};
                else
                    xtickLabelSet = {xtickLabelSet{:},tickLabel};
                end
            end
            xticks(xtickSet);
            xticklabels(xtickLabelSet);

            ytickSet = sort(unique(round(activeForceLengthCurveAnnotationPoints.y',1)));
            ytickLabelSet={''};
            for k=1:1:length(ytickSet)
                if(ytickSet(1,k)==0)
                    tickLabel = '0';
                else
                    tickLabel = sprintf('%1.1f%s',ytickSet(1,k),'$$f_o^M$$');
                end
                
                if(k==1)
                    ytickLabelSet={tickLabel};
                else
                    ytickLabelSet = {ytickLabelSet{:},tickLabel};
                end
            end
            yticks(ytickSet);
            yticklabels(ytickLabelSet);            


        else
            xtickSet = round(sort(unique(activeForceLengthCurveAnnotationPoints.x')),2);
            xtickSubSet= [xtickSet(1,1),xtickSet(1,2),1,...
                          round(lpeNOne,2),xtickSet(1,end)];
            xtickSubSet=round(sort(unique(xtickSubSet)),2);
            xticks(xtickSubSet);    
            xlabel('Norm. Length ($$\ell/\ell^{M}_o$$)');
            ylabel('Norm. Force ($$f/f^{M}_o$$)');            

        end

        legend('Location','NorthWest');
        legend boxoff;
        axis tight;
        ylim([0,1.05]);
        box off;

    subplot('Position',reshape(subPlotPanel(2,1,:),1,4));

       

      
        %Passive-force length plots
        fill([curveSampleECMHalf.x;...
              curveSampleECMHalf.x(end,1);...
              curveSampleECMHalf.x(1,1)].*2,...
             [curveSampleECMHalf.y;...
              curveSampleECMHalf.y(1,1);...
              curveSampleECMHalf.y(1,1)].*ecmForceFraction,...
              colorEcm,...
              'EdgeColor','none',...
              'FaceAlpha',0.5,...
              'HandleVisibility','off');
        hold on;


        %Active titin forces
        plot(curveSampleTitinActive.x.*2,...
             curveSampleTitinActive.y.*(1-ecmForceFraction) ...
             +curveSampleECMHalf.y.*(ecmForceFraction),...
             '-','Color',colorTitinActive,...
                'LineWidth',lineWidthModel,...
                'DisplayName',labelModel,...
                'HandleVisibility','off');
        hold on;            
        idxMid = 1;
        titinActive = ...
            curveSampleTitinActive.y.*(1-ecmForceFraction)...
           +curveSampleECMHalf.y.*(ecmForceFraction);
        while titinActive(idxMid,1) < 0.25
            idxMid=idxMid+1;
        end

        text(curveSampleTitinActive.x(idxMid,1).*2-0.025,...
             titinActive(idxMid,1).*(1-ecmForceFraction),...
             'Model: Active titin \& ECM',...
             'HorizontalAlignment','left',...
             'VerticalAlignment','bottom',...
             'FontSize',6,'Rotation',75);
        hold on;


        plot([lpeNOne;lpeNOne+0.05],[1;1],'k','HandleVisibility','off');
        hold on;
        plot([lpeNOne;lpeNOne+0.05],[1.01;1.01].*ecmForceFraction,...
            'k','HandleVisibility','off');
        hold on;
        plot([lpeNOne;lpeNOne+0.05],[0.99;0.99].*ecmForceFraction,...
            'k','HandleVisibility','off');
        hold on;

        plot([lpeNOne+0.05;lpeNOne+0.05],[1.01*ecmForceFraction;1],...
            'k','HandleVisibility','off');
        hold on;
        plot([lpeNOne+0.05;lpeNOne+0.05],[0;0.99*ecmForceFraction],...
            'k','HandleVisibility','off');
        hold on;

        plot([lpeNOne;lpeNOne+0.025],[1;1].*0,'k',...
            'HandleVisibility','off');
        hold on;
        text(lpeNOne+0.075, 0.5*ecmForceFraction,...
             'Model: ECM, $$f^{ECM}$$',...
             'FontSize',6,'Rotation',45);
        hold on;
        text(lpeNOne+0.075, ecmForceFraction + (1-ecmForceFraction).*0.5,...
             'Model: Passive Titin, $$f^{1},\,f^{2}$$',...
             'FontSize',6,'Rotation',45);
        hold on;

        titinTopX = [curveSampleTitin.x;...
                     fliplr(curveSampleECMHalf.x')'].*2;

        titinTopY = [(curveSampleTitin.y.*(1-ecmForceFraction)...
                      +curveSampleECMHalf.y.*ecmForceFraction);...
                     fliplr((curveSampleECMHalf.y.*ecmForceFraction)')'];

        fill(titinTopX,...
             titinTopY,...
             colorTitin,...
             'EdgeColor','none',...
             'FaceAlpha',0.5,...
              'HandleVisibility','off');
        hold on;

        plot(curveSampleTitin.x.*2 ,...
             (curveSampleTitin.y.*(1-ecmForceFraction) ...
             +curveSampleECMHalf.y.*ecmForceFraction),...
             '-','Color',colorPassive,...
                'LineWidth',lineWidthModel,...
                'DisplayName',[labelModel,'$$f^{PE}$$']);
        hold on

        yticks([0,ecmForceFraction,1]);
        yticklabels({'0','$$\lambda^{ECM} f^M_o$$','$$f^M_o$$'})

        xticks([lpeNZero,lpeNOne]);
        xticklabels({'\ell^M_{s}','\ell^M_{p}'});


        legend('Location','NorthWest');
        legend boxoff;
        axis tight;
        ylim([0,1.05]);
        box off;



    subplot('Position',reshape(subPlotPanel(1,2,:),1,4));

        if(flag_addCrossbridgeStrain==1)
            plot(fvCalSample.x,fvCalSample.y,...
                '-','Color',colorModel,...
                    'LineWidth',lineWidthModel,...
                    'DisplayName',[labelModel,' $$f^M_v$$']);
            hold on;

        else
            plot(fvSample.x,fvSample.y,...
                '-','Color',colorModel,...
                    'LineWidth',lineWidthModel,...
                    'DisplayName',[labelModel,' $$f^M_v$$']);
            hold on;
        end

        if(flag_addMat156==1)
            plot(fvSample.x,fvSample.y,'--','Color',colorMAT156,...
                'LineWidth',lineWidthMAT156,...
                'DisplayName',[labelMAT156,' $$f^M_v$$']);
            hold on;
        end

        plot(0,1,'.','Color',[0,0,0], 'HandleVisibility','off');
        hold on;
        plot([-1,0,0],[1,1,0],'--',...
            'Color',[1,1,1].*0.5,...
            'HandleVisibility','off');
        hold on;

        vmax = muscleModelStruct.musculotendon.maximumNormalizedFiberVelocity;
        vmaxStr = sprintf('%1.2f',vmax);
        vmaxHalfStr = sprintf('%1.2f',vmax*0.5);
        
        vmaxNegStr = ['-',vmaxStr,'$$\ell^{M}_o$$'];
        vmaxNegHalfStr = ['-',vmaxHalfStr,'$$\ell^{M}_o$$'];        
        vmaxPosStr = ['+',vmaxStr,'$$\ell^{M}_o$$'];
      
        if(flag_genericCurveConfig==1)
            vmaxNegStr = '$$-v^M_{max}$$';
            vmaxNegHalfStr = '$$-0.5 v^M_{max}$$';            
            vmaxPosStr = '$$v^M_{max}$$';
        end


        xticks([-1,-0.5,0,1]);
        xticklabels({vmaxNegStr,vmaxNegHalfStr,'0',vmaxPosStr});

        fvHalf = calcQuadraticBezierYFcnXDerivative(-0.5,...
            normMuscleQuadraticCurves.fiberForceVelocityCurve,...
            0);
        fvE     = muscleModelStruct.musculotendon.forceVelocityMultiplierAtLowEccentricFiberVelocity;
        fvEMax  = muscleModelStruct.musculotendon.forceVelocityMultiplierAtMaximumEccentricFiberVelocity;        

        if(flag_genericCurveConfig==1)
            ytickSet = [0,round(fvHalf,1),1.0,round(fvE,1),round(fvEMax,1)];
            ytickLabelSet={''};
            for k=1:1:length(ytickSet)
                if(ytickSet(1,k)==0)
                    tickLabel = '0';
                else
                    tickLabel = sprintf('%1.1f%s',ytickSet(1,k),'$$f_o^M$$');
                end
                
                if(k==1)
                    ytickLabelSet={tickLabel};
                else
                    ytickLabelSet = {ytickLabelSet{:},tickLabel};
                end
            end
            yticks(ytickSet);
            yticklabels(ytickLabelSet); 
        else
            yticks([0,round(fvHalf,2),1.0,round(fvE,2),round(fvEMax,2)]);
        end
        plot(-0.5,fvHalf,'.','Color',[0,0,0], 'HandleVisibility','off');
        hold on;
        plot([-1,-0.5,-0.5],[1,1,0].*fvHalf,'--','Color',[1,1,1].*0.5, 'HandleVisibility','off');
        hold on;        

        box off;
        axis tight;

        if(flag_genericCurveConfig==0)
            xlabel('Norm. Velocity ($$v/v^{M}_{max}$$)');
            ylabel('Norm. Force ($$f/f^{M}_o$$)');
        end
        title('B. Force-Velocity Relation');        
        
    subplot('Position',reshape(subPlotPanel(1,3,:),1,4));  
        plot(ftSample.x,ftSample.y,...
            '-','Color',colorModel,...
                'LineWidth',lineWidthModel,...
                'DisplayName',[labelModel,' $$f^T$$']);
        hold on;

        lt0 = 1;
        lt1 = normMuscleQuadraticCurves.tendonForceLengthCurve.xEnd(1,2);
        lt2 = lt0 + muscleModelStruct.musculotendon.tendonStrainAtOneNormForce;

        ft0 = 0;
        ft1 = normMuscleQuadraticCurves.tendonForceLengthCurve.yEnd(1,2);        
        ft2 = 1;

        
        plot(lt0,ft0,'.','Color',[0,0,0],'HandleVisibility','off');
        hold on;
        plot(lt1,ft1,'.','Color',[0,0,0],'HandleVisibility','off');
        hold on;
        plot(lt2,ft2,'.','Color',[0,0,0],'HandleVisibility','off');
        hold on;

        if(flag_genericCurveConfig==1)
            xtickSet = sort(unique(round([lt0,lt1,lt2],3)));
            xtickLabelSet={'$$\ell^T_s$$','$$\ell^T_{toe}$$','$$\ell^T_{o}$$'};
            xticks(xtickSet);
            xticklabels(xtickLabelSet);

            ytickSet = sort(unique(round([ft0,ft1,ft2]',3)));
            ytickLabelSet={'$$0$$','$$f^T_{toe}$$','$$f^M_{o}$$'};
            yticks(ytickSet);
            yticklabels(ytickLabelSet);            


        else
            xticks(round([lt0,lt1,lt2],3));
            yticks(round([ft0,ft1,ft2],3));          
            xlabel('Norm. Length ($$\ell/\ell^{T}_{s}$$)');
            ylabel('Norm. Force ($$f/f^{M}_o$$)');
        end

        hold on;
        box off;

        title('C. Tendon Force-Length Relation');
        here=1;

        axis tight;
        ylim([0,1.05]);

