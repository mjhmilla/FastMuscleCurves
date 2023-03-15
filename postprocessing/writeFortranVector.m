function success = writeFortranVector(dataX, dataY, tableId, fullFilePath)

success = 0;

assert(size(dataX,2)==1 && size(dataY,2)==1,'Error: both dataX and dataY must be column vectors');
assert(size(dataX,1)==size(dataY,1),'Error: dataX and dataY must have the same length');
assert(abs(round(tableId)-tableId) < eps,'Error tableId must be an integer');

fid = fopen(fullFilePath,'w');
assert(fid,'Error: fullFilePath failed to open');

idStr = sprintf('%i',tableId);
for i=length(idStr):1:9
	idStr = [' ',idStr];
end	


fprintf(fid, '$\n');
fprintf(fid,'*DEFINE_CURVE\n');
fprintf(fid,'$#    lcid      sidr       sfa       sfo      offa      offo    dattyp\n');
fprintf(fid,'%s         0       1.0       1.0       0.0       0.0\n',idStr);
fprintf(fid,'$#                a1                  o1\n');

maxOrderX =round(log10(max(abs(dataX))));
maxOrderY =round(log10(max(abs(dataY))));

orderMax = max([maxOrderX,maxOrderY]);

fieldWidth = 20;%As in 20 characters

if(orderMax < 1)
    orderMax=1;
end

printPattern = ['%1.',num2str(fieldWidth-4-orderMax),'f'];


for i=1:1:size(dataX,1)

	dataXStr = sprintf(printPattern,dataX(i,1));
	dataYStr = sprintf(printPattern,dataY(i,1));

	for j=length(dataXStr):1:(fieldWidth-1)
		dataXStr = [' ',dataXStr];
    end
	for j=length(dataYStr):1:(fieldWidth-1)
		dataYStr = [' ',dataYStr];
    end

    fprintf(fid,'%s%s\n',dataXStr,dataYStr);

end




fclose(fid);

success=1;