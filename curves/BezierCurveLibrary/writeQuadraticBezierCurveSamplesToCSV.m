function status = writeQuadraticBezierCurveSamplesToCSV(...
    structOfMuscleCurves,nPts,muscleName,folder)

status=0;

curveNames = fieldnames(structOfMuscleCurves);

for indexCurve=1:1:length(curveNames)
    fid = fopen([folder,muscleName,'_',curveNames{indexCurve},'.csv'],'w');

    fprintf(fid,'arg,val,der1,der2\n');

    xEnd = structOfMuscleCurves.(curveNames{indexCurve}).xEnd;

    xSpan = diff(xEnd);
    x0 = xEnd(1,1) - xSpan*0.1;
    x1 = xEnd(1,2) + xSpan*0.1;

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