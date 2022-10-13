function success = writeBezierCurveStructuresToFortran(structOfMuscleCurves,...
    muscleName,folder)

success=0;
curveNames = fieldnames(structOfMuscleCurves);

for indexCurve=1:1:length(curveNames)
    if(contains(curveNames{indexCurve},'use')==0)
        fid = fopen([folder,muscleName,'_',curveNames{indexCurve},'.f'],'w');
    
        xpts = structOfMuscleCurves.(curveNames{indexCurve}).xpts;
        ypts = structOfMuscleCurves.(curveNames{indexCurve}).ypts;
        
        xEnd = structOfMuscleCurves.(curveNames{indexCurve}).xEnd;
        yEnd = structOfMuscleCurves.(curveNames{indexCurve}).yEnd;
    
        dydxEnd = structOfMuscleCurves.(curveNames{indexCurve}).dydxEnd;
        d2ydx2End = structOfMuscleCurves.(curveNames{indexCurve}).d2ydx2End;
        
    %     fprintf(fid,'xpts\n');
    %     for i=1:1:size(xpts,1)
    %         for j=1:1:size(xpts,2)
    %             if(j < size(xpts,2))
    %                 fprintf(fid,['%1.16e',delimiter],xpts(i,j));
    %             else
    %                 fprintf(fid,['%1.16e',delimiter],xpts(i,j));
    %             end
    %         end
    %     end
    
        fprintf(fid,'      real*8, dimension(%d,%d)::xPts \n',   size(xpts,1),size(xpts,2));
        fprintf(fid,'      real*8, dimension(%d,%d)::yPts \n',   size(ypts,1),size(ypts,2));
        fprintf(fid,'      real*8, dimension(%d)::xEnd \n',      size(xEnd,2));
        fprintf(fid,'      real*8, dimension(%d)::yEnd \n',      size(yEnd,2));
        fprintf(fid,'      real*8, dimension(%d)::dydxEnd \n',   size(dydxEnd,2));
        fprintf(fid,'      real*8, dimension(%d)::d2ydx2End \n', size(d2ydx2End,2));
        fprintf(fid,'      ncol = %d\n',size(xpts,2));

        fid=writeArrayToFortran(fid,xpts,'xPts');
        fid=writeArrayToFortran(fid,ypts,'yPts');
        fid=writeArrayToFortran(fid,xEnd,'xEnd');
        fid=writeArrayToFortran(fid,yEnd,'yEnd');
        fid=writeArrayToFortran(fid,dydxEnd,    'dydxEnd');
        fid=writeArrayToFortran(fid,d2ydx2End,  'd2ydx2End');
        fclose(fid);
    end
end

success=1;