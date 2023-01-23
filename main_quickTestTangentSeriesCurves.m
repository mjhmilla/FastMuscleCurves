clc;
close all;
clear all;

x = [0:0.01:1]';

dy  = tanh(x);
y   = log(cosh(x));

x0 = 0;
c0 = 0;
c1 = log(cosh(x0));
c2 = tanh(x0);

sechxSq = sech(x0)*sech(x0);
tanhxSq = tanh(x0)*tanh(x0);

c3 = sechxSq;
c4 = -2*sechxSq*tanh(x0);
c5 = -2*( -2*sechxSq*tanhxSq + sechxSq*sechxSq);

iy  = c0 ...
        + (c1).*x ...
        + (1/2).*(c2).*(x.^2)...
        + (1/3).*(c3).*(x.^3)...
        + (1/4).*(c4).*(x.^4)...
        + (1/5).*(c5).*(x.^5);

dyfcn = @(argT,argY)tanh(argT);
yfcn = @(argT,argY)log(cosh(argT));



[tyNum, yNum] = ode45(dyfcn,x,y(1,1));
[tiyNum, iyNum] = ode45(yfcn,x,iy(1,1));

%yNum=yNumSol.y;
%iyNum=iyNumSol.y;


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
    plot(x,abs(iy),'-','Color',[0,0,0]);
    hold on;
    plot(x,abs(iyNum),'--','Color',[0,0,1]);
    hold on;

    box off;

    xlabel('X');
    ylabel('Y');
    title('int( int( tanh(x) ))');

    subplot(2,3,6);
        plot(x,iy-iyNum,'-','Color',[1,0,0]);
        hold on;
    
        box off;
    
        xlabel('X');
        ylabel('Y');
        title('Err: int( int( tanh(x) ))');