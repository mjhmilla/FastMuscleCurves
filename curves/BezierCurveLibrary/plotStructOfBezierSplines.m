function structOfFigures = plotStructOfBezierSplines( structOfFigures, ...
                                         curveStruct ,...
                                         ignoreStructsWithThisKeyWord,...
                                         colorOverride, lineWidthOverride)


curveNames =fieldnames(curveStruct);

%structOfFigures  = [];
structOfFiguresFields = [];
if(isempty(structOfFigures)==0)
    structOfFiguresFields = fieldnames(structOfFigures);
end

for i=1:1:length(curveNames)
  
  idxKeyWord = strfind(curveNames{i},ignoreStructsWithThisKeyWord);
  
  if(isempty(curveStruct.(curveNames{i})) == 0 ...
     && isempty(idxKeyWord) == 1)
    
    if isempty(structOfFiguresFields)
        structOfFigures.(curveNames{i}) = figure;
    elseif isNameInList(curveNames{i}, structOfFiguresFields)==0
        structOfFigures.(curveNames{i}) = figure;
    else
        figure(structOfFigures.(curveNames{i}));
    end

    if(strcmp(curveNames{i},'fiberForceVelocityCurve')==1)
        here=1;
    end

    curveSample = calcBezierYFcnXCurveSampleVector(...
                    curveStruct.(curveNames{i}), 500,[]);

    xmin = min(curveSample.x);
    xmax = max(curveSample.x);
    ymin = min(curveSample.y);
    ymax = max(curveSample.y);

    xV   = curveSample.x;
    yV   = curveSample.y;
    y1V  = curveSample.dydx;
    y2V  = curveSample.d2ydx2;

    subplot(2,2,1);

    yColor = [0,0,0];
    
    if(isempty(colorOverride)==0)
        yColor = colorOverride;        
    end
    yLineWidth=1;
    if(isempty(lineWidthOverride)==0)
        yLineWidth = lineWidthOverride;        
    end

    plot(curveSample.x, curveSample.y,'Color',yColor,...
        'LineWidth',yLineWidth);
      hold on;        

      xlabel('x');
      ylabel('y');
      title(curveNames{i});
      grid on;
      axis square;
      box off;            
      xlim([xmin,xmax]);

    subplot(2,2,2);   

    dydxColor = [1,0,0];
    if(isempty(colorOverride)==0)
        dydxColor = colorOverride;
    end

    dydxLineWidth=1;
    if(isempty(lineWidthOverride)==0)
        dydxLineWidth = lineWidthOverride;        
    end
    
    plot(curveSample.x, curveSample.dydx,'Color',dydxColor,...
        'LineWidth',dydxLineWidth);
      hold on;            
      xlabel('x');
      ylabel('dy/dx');
      grid on;
      axis square;
      box off;            
      xlim([xmin,xmax]); 

    subplot(2,2,3);   

    d2ydx2Color = [0,0,1];
    if(isempty(colorOverride)==0)
        d2ydx2Color = colorOverride;
    end
    
    d2ydx2LineWidth=1;
    if(isempty(lineWidthOverride)==0)
        d2ydx2LineWidth = lineWidthOverride;        
    end

    plot(curveSample.x, curveSample.d2ydx2,'Color',d2ydx2Color,...
         'LineWidth',d2ydx2LineWidth);
      hold on;            
      xlabel('x');
      ylabel('d2y/dx2');
      grid on;
      axis square;
      box off;            
      xlim([xmin,xmax]); 

    subplot(2,2,4);   

    d3ydx3Color = [0,1,0];
    if(isempty(colorOverride)==0)
        d3ydx3Color = colorOverride;
    end

    d3ydx3LineWidth=1;
    if(isempty(lineWidthOverride)==0)
        d3ydx3LineWidth = lineWidthOverride;        
    end

    plot(curveSample.x, curveSample.d3ydx3,'Color',d3ydx3Color,...
         'LineWidth',d3ydx3LineWidth);
      hold on;            
      xlabel('x');
      ylabel('d3y/dx2');
      grid on;
      axis square;
      box off;            
      xlim([xmin,xmax]);       

%     if(isempty(curveStruct.(curveNames{i}).integral)==0)  
%         intYdx = curveSample.intYdx;
% 
%         subplot(2,2,4)
%             plot(curveSample.x, intYdx,'g');
%             hold on;
%             xlabel('x');
%             ylabel(['int(y)']);
%             axis square;
%             box off;                
%             grid on;
%             hold on;
%             xlim([xmin,xmax]);  
% 
%     end
  end
  
end


