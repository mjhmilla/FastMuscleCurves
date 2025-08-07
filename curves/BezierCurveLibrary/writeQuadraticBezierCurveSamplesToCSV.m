function status = writeQuadraticBezierCurveSamplesToCSV(...
    structOfMuscleCurves,nPts,muscleName,folder,flag_extendedDomain)

status=0;

curveNames = fieldnames(structOfMuscleCurves);

for indexCurve=1:1:length(curveNames)

    if(length(muscleName)>0)
        fid = fopen([folder,muscleName,'_',curveNames{indexCurve},'.csv'],'w');
    else
        fid = fopen([folder,curveNames{indexCurve},'.csv'],'w');
    end

    fprintf(fid,'arg,val,der1,der2\n');

    xEnd = structOfMuscleCurves.(curveNames{indexCurve}).xEnd;

    xSpan = diff(xEnd);
    x0 = xEnd(1,1) - xSpan*0.1;
    x1 = xEnd(1,2) + xSpan*0.1;

    if(flag_extendedDomain==1)
        dx = x1-x0;
        x0 = x0-3*dx;
        x1 = x1+3*dx;
        nPts=nPts*3;
    end

    xSamples = [x0:((x1-x0)/(nPts)):x1]';

    for i=1:1:length(xSamples)

        x       = xSamples(i,1);
        y       = calcQuadraticBezierYFcnXDerivative(x, ...
                    structOfMuscleCurves.(curveNames{indexCurve}), 0);
        dydx    = calcQuadraticBezierYFcnXDerivative(x, ...
                    structOfMuscleCurves.(curveNames{indexCurve}), 1);
        d2ydx2  = calcQuadraticBezierYFcnXDerivative(x, ...
                    structOfMuscleCurves.(curveNames{indexCurve}), 2);
        fprintf(fid,'%e,%e,%e,%e\n',x,y,dydx,d2ydx2);        

    end
    fclose(fid);
end



status=1;