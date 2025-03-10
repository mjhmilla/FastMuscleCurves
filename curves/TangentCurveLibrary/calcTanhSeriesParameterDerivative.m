function partialDerivative = calcTanhSeriesParameterDerivative(x,coeffs,...
    functionDerivativeOrder, indexOfParameterDerivative)

assert(functionDerivativeOrder==0 || functionDerivativeOrder==1,...
       'Error: this function is valid only for 0 or 1st order derivatives');

assert(size(coeffs,1) == 1,...
    ['Error: this function only evaluates the',...
    ' partial derivative for a single segment']);

validIndex=0;
for i=1:1:6
    if(indexOfParameterDerivative==i)
        validIndex=1;
    end
end

assert(validIndex==1,'Error: indexOfParameterDerivative must be 1-6');

partialDerivative = 0;


%derIndex = find(derivativeOrder(1,:)>0);

A = coeffs(1,1);
B = coeffs(1,2);
C = coeffs(1,3);
D = coeffs(1,4);
E = coeffs(1,5);
F = coeffs(1,6);

if(functionDerivativeOrder==0)
    switch indexOfParameterDerivative
        case 1
    %         d/dA
    %         arg = (x-B)/C;
    %         logCoshArg = log(1.+exp(-2*arg))-log(2)+arg;    
    %         y = y + (D*(x-B) + A*C*(logCoshArg) + C + E);
    %
    %         Darg_DA = 0
    %         DlogCoshArg_DA = 0
    %         Dy_DA = C*(logCoshArg)
            
            arg = (x-B)/C;

            argMin = -log(realmax)*0.1;
            logCoshArg = 0;
            if(arg <= argMin)
                logCoshArg = -arg-log(2);
            else
                logCoshArg = log(1.+exp(-2*arg))-log(2)+arg;
            end            
            %y = y + (D*(x-B) + A*C*(logCoshArg) + C + E);
    
            partialDerivative = C*logCoshArg;
        case 2
    %         d/dB
    %         arg = (x-B)/C;
    %         logCoshArg = log(1.+exp(-2*arg))-log(2)+arg;    
    %         y = y + (D*(x-B) + A*C*(logCoshArg) + C + E);
    %
    %         Darg_DB = -1/C
    %         DlogCoshArg_DB = 0
    %         Dy_DB = C*(logCoshArg)
              arg = (x-B)/C;
              argMin = -log(realmax)*0.1;

              Darg_DB = -1/C;
	          DlogCoshArg_DB =0;               

              if(arg <= argMin)
                  DlogCoshArg_DB = (2/C) +Darg_DB; 
              else
                  DlogCoshArg_DB =( 2*exp(-2*arg) / (C*(exp(-2*arg)+1)) ) +Darg_DB; 
              end 

	                   
              partialDerivative = -D + A*C*DlogCoshArg_DB;
        case 3
    %         d/dC
    %         arg = (x-B)/C;
    %         logCoshArg = log(1.+exp(-2*arg))-log(2)+arg;    
    %         y = y + (D*(x-B) + A*C*(logCoshArg) + C + E);
    %
    %         Darg_DC = -(x-B)/(C*C)
    %         DlogCoshArg_DC = (-2*exp(-2*arg)*Darg_DC) / (1.+exp(-2*arg)) + Darg_DC
    %         Dy_DC = A*logCoshArg + A*C*DlogCoshArg_DC + 1   
            arg = (x-B)/C;
            argMin = -log(realmax)*0.1;
            logCoshArg = 0;
            DlogCoshArg_DC = 0;
            Darg_DC = -(x-B)/(C*C);            
            if(arg <= argMin)
                logCoshArg = -arg-log(2);
                DlogCoshArg_DC = (-2*Darg_DC) + Darg_DC;                
            else
                logCoshArg = log(1.+exp(-2*arg))-log(2)+arg;  
                DlogCoshArg_DC = (-2*exp(-2*arg)*Darg_DC) / (1.+exp(-2*arg)) + Darg_DC;
            end


            partialDerivative = A*logCoshArg + A*C*DlogCoshArg_DC + 1;
    
        case 4
            %d/dD
    %         arg = (x-B)/C;
    %         logCoshArg = log(1.+exp(-2*arg))-log(2)+arg;    
    %         y = y + (D*(x-B) + A*C*(logCoshArg) + C + E);
    %
    %         Darg_DD = 0
    %         DlogCoshArg_DD = 0
    %         Dy_DD = x-B       
            partialDerivative = x-B;
        case 5
            %d/dE
    %         arg = (x-B)/C;
    %         logCoshArg = log(1.+exp(-2*arg))-log(2)+arg;    
    %         y = y + (D*(x-B) + A*C*(logCoshArg) + C + E);
    %
    %         Darg_DE = 0
    %         DlogCoshArg_DE = 0
    %         Dy_DE = 1          
              partialDerivative=1;
    
        case 6
            %d/dF
    %         arg = (x-B)/C;
    %         logCoshArg = log(1.+exp(-2*arg))-log(2)+arg;    
    %         y = y + (D*(x-B) + A*C*(logCoshArg) + C + E);
    %
    %         Darg_DF = 0
    %         DlogCoshArg_DF = 0
    %         Dy_DF = 0  
              partialDerivative=0;
    
        otherwise
            assert(0,'Error: partial derivative index must be [1-6]');
    end
end

if(functionDerivativeOrder==1)
    switch indexOfParameterDerivative
        case 1
    %         d/dA
    %         dy = A*tanh((x-B)/C) + D;
    %         Ddy_DA = tanh((x-B)/C)
    
            partialDerivative = tanh((x-B)/C);

        case 2
    %         d/dB
    %         dy = (A*tanh((x-B)/C) + D);
    %         

            sechArg = sech((x-B)/C);
            partialDerivative = -A*(sechArg*sechArg)/C;

        case 3
    %         d/dC
    %         dy = (A*tanh((x-B)/C) + D);              
            sechArg = sech((x-B)/C);
            partialDerivative = -A*(x-B)*sechArg*sechArg / (C*C);
    
        case 4
    %         d/dD
    %         dy = (A*tanh((x-B)/C) + D);
   
            partialDerivative = 1;
        case 5
    %         d/dE
    %         dy = (A*tanh((x-B)/C) + D);
         
              partialDerivative=0;
    
        case 6
    %         d/dF
    %         dy = (A*tanh((x-B)/C) + D);

              partialDerivative=0;
    
        otherwise
            assert(0,'Error: partial derivative index must be [1-6]');
    end
end


