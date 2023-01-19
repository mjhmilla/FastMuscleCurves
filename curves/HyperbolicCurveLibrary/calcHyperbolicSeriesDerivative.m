function value = calcHyperbolicSeriesDerivative(x, coeffs, derivativeOrder)

value = nan;

switch derivativeOrder
    case -1
        inty = 0;

        for i=1:1:size(coeffs,1)
            A = coeffs(i,1);
            B = coeffs(i,2);
            C = coeffs(i,3);
            D = coeffs(i,4);
            E = coeffs(i,5);
            F = coeffs(i,6);
            x2=x*x;
            arg = abs((x-B)/C);

            x2=x*x;
            arg = abs((x-B)/C);
            
            int_lncoshx=arg*arg*0.5-log(2)*arg; %Limit, apparently
%            int_lncoshx=arg*arg*0.5-log(2)*arg;
%            for j=1:1:10
%                int_lncoshx = int_lncoshx + ((-1)^(j))*(exp(-2*j*arg)/(2*j*j));
%            end

           inty = inty + ((C+E-B*D)*x ...
                       + D*x2*0.5 ...
                       + A*C*(int_lncoshx)...
                       + F);            
         end
 
         value=inty;
    case 0
        y = 0;
        for i=1:1:size(coeffs,1)
            A = coeffs(i,1);
            B = coeffs(i,2);
            C = coeffs(i,3);
            D = coeffs(i,4);
            E = coeffs(i,5);


            y = y + (-B*D + D*x + A*C*log( cosh((x-B)/C) ) + C + E);
        end
        value=y;
    case 1
        dy = 0;
        for i=1:1:size(coeffs,1)
            A = coeffs(i,1);
            B = coeffs(i,2);
            C = coeffs(i,3);
            D = coeffs(i,4);

            dy = dy + (A*tanh((x-B)/C) + D);
        end
        value=dy;
        
    otherwise
        assert(0,'Error: derivative order must be [-1,0,1]');
end