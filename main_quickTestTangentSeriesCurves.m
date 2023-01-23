clc;
close all;
clear all;

x = [0:0.01:1]';

dy  = tanh(x)+1;
y   = log(cosh(x))+x;
iy  = 2.*atan(exp(x)) ...
     + 0.5.*(2.*x.*log(exp(2.*x)+1) - 3.*x.*x - 2.*x.*log(2))  ...
     + x.*x.*0.5;


yNum = cumtrapz(x,dy,1) + y(1,1);
iyNum= cumtrapz(x,y,1) + iy(1,1);



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