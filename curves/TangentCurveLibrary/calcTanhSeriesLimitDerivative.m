function value = calcTanhSeriesLimitDerivative(x, coeffs, derivativeOrder)

value = nan;

switch derivativeOrder

    case 0
        y = 0;
        for i=1:1:size(coeffs,1)
            A = coeffs(i,1);
            B = coeffs(i,2);
            C = coeffs(i,3);
            D = coeffs(i,4);
            E = coeffs(i,5);
            
            %logCoshArg = log(1.+exp(-2*arg))-log(2)+arg;
            %
            % As arg -> +inf
            %
            % log(1.+exp(-2*arg)) -> log(exp(-2*arg)) -> -2*arg
            %
            % Thus
            %
            % logCoshArg -> -2*arg + arg - log(2)
            %            -> -arg
            %
            %And so
            %
            % y = y + (D*(x-B) + A*C*(logCoshArg) + C + E)
            % 
            % since
            %
            % arg = (x-B)/C
            %
            %y = y + (D*C*arg - A*C*(arg) + C + E)
            %
            %y = y + (C*arg*(D-A) + C + E)
    
            arg = (x-B)/C;            
            y = y + C*arg*(D-A) + C + E;
        end
        value=y;
    case 1
        dy = 0;
        for i=1:1:size(coeffs,1)
            A = coeffs(i,1);
            B = coeffs(i,2);
            C = coeffs(i,3);
            D = coeffs(i,4);    

            % This function can be evaluated for large values of (x-B)/C
            %
            %
            %

            dy = dy + (A*tanh((x-B)/C) + D);
        end
        value=dy;
        
    otherwise
        assert(0,'Error: derivative order must be [0,1]');
end

value=real(value);