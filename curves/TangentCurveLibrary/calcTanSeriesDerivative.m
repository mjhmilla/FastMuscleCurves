function value = calcTanSeriesDerivative(x, coeffs, derivativeOrder)

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
     
            argXBC = (x-B)/C;
            
            dC = 1/C;
            C2 = C*C;
            C3 = C2*C;
            C4 = C3*C;
            B2 = B*B;
            
            x2 = x*x;
            
            atanXBC = atan(argXBC);
            
            inty = inty +...
                A*C*(...
                  ((x2*0.5-B*x) * atanXBC ...
                      -(...
                          (...
                            (-C4-B2*C2)*atanXBC...
                          )*(0.5*dC) ...
                          +(C2*x)*0.5 ...
                        )*dC ...
                  )*dC ...
                 -( ...
                    x*log(argXBC*argXBC+1) ...
                    -(2*(...
                          (B*C2*log(x2-2*B*x+C2+B2))*0.5 ...
                          -C3*atanXBC+C2*x ...
                        ) ...
                      )*(dC*dC) ...
                  )*0.5 ...
                )...
                +(D*x2)*0.5...
                +E*x ...
                +F;            
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
            argXB = x-B;
            argXBC = argXB/C;
            
            y = y+ ...
                A*( argXB*atan(argXBC)-log(argXBC*argXBC+1)*0.5*C)+D*x + E;

        end
        value=y;
    case 1
        dy = 0;
        for i=1:1:size(coeffs,1)
            A = coeffs(i,1);
            B = coeffs(i,2);
            C = coeffs(i,3);
            D = coeffs(i,4);            
            dy = dy + A*atan((x-B)/C) + D;
        end
        value=dy;
        
    otherwise
        assert(0,'Error: derivative order must be [-1,0,1]');
end