function value = calcTanhSeriesDerivative(x, coeffs, derivativeOrder)

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
            C2 = C*C;
            
            inty=inty+...
                0.5.*((A*C2).*polylog(2,-exp((2.*(B-x))./C)) ...
                     + A.*(B-x).*(B+C*log(4)-x) ...
                     + x.*(2.*(-B*D+C+E) + D.*x)...
                     ) ...
               + F;         
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
            arg = (x-B)/C;
            logCoshArg = log(1.+exp(-2*arg))-log(2)+arg;
    
            y = y + (D*(x-B) + A*C*(logCoshArg) + C + E);
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

value=real(value);