clc;
close all;
clear all;

x = [-2:0.01:2]';


dy  = tanh(x);
y   = log(cosh(x));
iy = 0.5.*(polylog(2,-exp(-2.*x)) ...
          - x.*(x + 2.*log(exp(-2.*x) + 1) ...
          - 2.*log(cosh(x))));

dyfcn = @(argT,argY)tanh(argT);
yfcn = @(argT,argY)(log(cosh(argT)));



[tyNum, yNum]   = ode45(dyfcn,x,y(1,1));
[tiyNum, iyNum] = ode45(yfcn,x,(iy(1,1)));

assert(max((tyNum-x))<eps);
assert(max((tiyNum-x))<eps);

figPlot=figure;
subplot(2,3,1);
    plot(x,dy,'-','Color',[0,0,0]);

    box off;

    xlabel('X');
    ylabel('Y');
    title('tanh(x)');


subplot(2,3,2);
    plot(x,y,'-','Color',[0,0,0]);
    hold on;
    plot(x,yNum,'--','Color',[0,0,1]);
    hold on;

    box off;

    xlabel('X');
    ylabel('Y');
    title('int( tanh(x) )');

    subplot(2,3,5);
        plot(x,y-yNum,'-','Color',[1,0,0]);
        hold on;
        box off;
    
        xlabel('X');
        ylabel('Y');
        title('Err: int( tanh(x) )');

subplot(2,3,3);
    plot(x,(iy),'-','Color',[0,0,0]);
    hold on;
    plot(x,iyNum,'--','Color',[0,0,1]);
    hold on;

    box off;

    xlabel('X');
    ylabel('Y');
    title('int( int( tanh(x) ))');

    subplot(2,3,6);
        plot(x,(iy)-iyNum,'-','Color',[1,0,0]);
        hold on;
    
        box off;
    
        xlabel('X');
        ylabel('Y');
        title('Err: int( int( tanh(x) ))');