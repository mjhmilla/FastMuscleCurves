function val = calcQuadraticBezierYFcnXDerivative(x, curveParams, der)
%%
% This function takes the control points for the Bezier curves
% x(t) and f(t) and treats it as a function:
%
%        val = y(x)
%
% To do this we assume that
% 1. x(t) is a monotonic increasing function
% 2. d/dt x(t) > 0
%
% In addition, to maintain continuity we assume that the second derivative
% of the curve goes to 0 at the end points so that we can linearly
% extrapolate the curve outside of [xmin, xmax]
%
% @param x: value to 
%
% @param curveParams : a structure with the fields
%           .xpts
%           .ypts
%           .integral
%
%   xpts: n x m matrix of Bezier control points defining x(t). 
%
%              n x m matrix 
%              n: Bezier curve order + 1 
%              m: # spline sections
%                
%              Each column defines the control points used for a Bezier 
%              section. Thus a quintic Bezier spline with 3 sections will 
%              have a 6 x 3 matrix for xpts
%              
%   ypts: n x m matrix of Bezier control points defining y(t)
%
%   integral: if the integral curve has not been numerically computed, 
%             then this field will be empty. Otherwise it will have fields
%             of
%
%                      .xptsN : points the integral is evaluated at  
%                      .yptsN : numerical value of the integral
%                      .y1ptsN: '' integral's first derivative.
%                      .y2ptsN: '' integral's 2nd derivative.
%
% @param der: integer defining the order of the desired derivative. 
%             The value of der is restricted to:
%             
%             0: y = f(x)
%             1: dy/dx  
%             2: d^2y/dx^2 
%             3: d^3y/dx^3
%
% @returns y : d^n/dx^n y(x), the n^th derivative of the function y(x)
%%

val = NaN;

assert((der >= -1 && der <= 3),'der must be within [0,3]');

%xpts = curveParams.xpts;
%ypts = curveParams.ypts;

nrow = size(curveParams.xpts,1);
ncol = size(curveParams.xpts,2);
xmin = curveParams.xEnd(1,1);
xmax = curveParams.xEnd(1,2);
npts = size(cuveParams.xpts,1);

assert(npts == 3);

if( (x < xmin) || (x > xmax) )
    idxEnd = NaN;
    if(x <= xmin)
       idxEnd = 1; 
    else
       idxEnd = 2;
    end
    
    %evaluate the desired derivative of y at x
    switch der           
        case 0
            x0   = curveParams.xEnd(idxEnd);
            y0   = curveParams.yEnd(idxEnd);
            dydx = curveParams.dydxEnd(idxEnd);
                        
            val = dydx*(x-x0) + y0;             
        case 1
            dydx = curveParams.dydxEnd(idxEnd);
            val  = dydx;
        otherwise
            val = 0;
    end        
else
    %Find the spline section that x is in
    options.moduloRange = [];
    options.tol    = eps; 
    xInterval      = [curveParams.xpts(   1,:); ...
                      curveParams.xpts(nrow,:)];
    col =  calcIndex(x, xInterval, options);

    u = NaN;

    x0 = curveParams.xpts(1,col);
    x1 = curveParams.xpts(2,col);
    x2 = curveParams.xpts(3,col);

    % x    = (1-u)^2 x0 + 2*u(1-u)*x1 + u^2 x2
    %      = (1-2*u  +u^2)*x0
    %        (  2*u-2*u^2)*x1
    %        (        u^2)*x2
    % 0   = (x0-2*x1+x2)*u^2 + (-2*x0+2*x1)*u + (x0 - x)

    a = x0-2*x1+x2;
    b = -2*x0 +2*x1;
    c = x0-x;

    t0 = sqrt(b*b-4*a*c);
    u = (-b + t0)/(2*a);        
    if( u < 0 || u > 1)
        u = (-b - t0)/(2*a);
    end
    assert( u >=0 && u <= 1);
    
end

    %Evaluate the desired derivative of y   
    switch der
        case -1
             assert(0);
           
        case 0            
            u2 = u*u;
            v  = (1-u);
            v2 = v*v;

            y0 = curveParams.ypts(1,col);
            y1 = curveParams.ypts(2,col);
            y2 = curveParams.ypts(3,col);
            
            val = v2*y0 + 2*u*v*y1 + u2*y2;                    
        case 1

            du =1;

            u2 = u*u;
            du2= 2*u;
            
            v  = (1-u);
            dv = -1;
            
            v2 = v*v;
            dv2=2*v;

            x0 = curveParams.xpts(1,col);
            x1 = curveParams.xpts(2,col);
            x2 = curveParams.xpts(3,col);

            y0 = curveParams.ypts(1,col);
            y1 = curveParams.ypts(2,col);
            y2 = curveParams.ypts(3,col);
            
            %y = v2*y0 + 2*u*v*y1 + u2*y2;
            dy= dv2*y0 + (2*du*v+2*u*dv)*y1 + du2*y2;

            %x = v2*x0 + 2*u*v*x1 + u2*x2;
            dx= dv2*x0 + (2*du*v+2*u*dv)*x1 + du2*x2;

            val = dy/dx;
            
        case 2

            du=1;

            u2 = u*u;
            du2= 2*u;
            d2u2= 2;

            v  = (1-u);
            dv = -1;
            %d2v= 0;
            
            v2 = v*v;
            dv2=2*v;
            d2v2=2;

            x0 = curveParams.xpts(1,col);
            x1 = curveParams.xpts(2,col);
            x2 = curveParams.xpts(3,col);

            y0 = curveParams.ypts(1,col);
            y1 = curveParams.ypts(2,col);
            y2 = curveParams.ypts(3,col);
            
            %y = v2*y0 + 2*u*v*y1 + u2*y2;
            dy =  dv2*y0  + (2*du*v+2*u*dv)*y1 + du2*y2;
            d2y= d2v2*y0 + (4*du*dv)*y1 + d2u2*y2;

            %x = v2*x0 + 2*u*v*x1 + u2*x2;
            dx =  dv2*x0 + (2*du*v+2*u*dv)*x1 + du2*x2;
            d2x= d2v2*x0 + (4*du*dv)*x1 + d2u2*x2;
  
            val  = (d2y*dx - dy*d2x)/(dx*dx);

        case 3
            du=1;

            u2 = u*u;
            du2= 2*u;
            d2u2= 2;

            v  = (1-u);
            dv = -1;
            %d2v= 0;
            
            v2 = v*v;
            dv2=2*v;
            d2v2=2;

            x0 = curveParams.xpts(1,col);
            x1 = curveParams.xpts(2,col);
            x2 = curveParams.xpts(3,col);

            y0 = curveParams.ypts(1,col);
            y1 = curveParams.ypts(2,col);
            y2 = curveParams.ypts(3,col);
            
            %y =   v2*y0            + 2*u*v*y1   + u2*y2;
            dy =  dv2*y0  + (2*du*v+2*u*dv)*y1  + du2*y2;
            d2y= d2v2*y0        + (4*du*dv)*y1 + d2u2*y2;
            %d3y= 0;

            %x =   v2*x0 +           2*u*v*x1 +   u2*x2;
            dx =  dv2*x0 + (2*du*v+2*u*dv)*x1 +  du2*x2;
            d2x= d2v2*x0 +       (4*du*dv)*x1 + d2u2*x2;
            %d3x = 0;

            val  = -2*(d2y*dx - dy*d2x)*(2*d2x*dx)/(dx*dx*dx);
           
    end
            
    
end
    