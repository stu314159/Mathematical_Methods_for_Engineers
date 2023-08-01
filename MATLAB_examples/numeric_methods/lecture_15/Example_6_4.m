%% Example 6-4

clear
clc
close 'all'

%% data
x = [1 2 4 5 7];
y = [52 5 -5 -40 10];

% get the generated Lagrange Polynomial
F = genLagrangePolyInterp(x,y);

% test at f(3)

fprintf('F(3) = %g \n',F(3));
% plot the data and the interpolant

xMin = 0;
xMax = 8;
nX = 100;
xSpace = linspace(xMin,xMax,nX);

figure(1)
plot(xSpace,F(xSpace),...
    x,y,'sk','linewidth',2,'markersize',10);
grid on
title('Lagrange Interpolation','fontsize',18,'fontweight','bold');
xlabel('X','fontsize',16,'fontweight','bold');
ylabel('Y','fontsize',16,'fontweight','bold');
set(gca,'fontsize',12,'fontweight','bold');

%% case 2: "Witch of Agnesi" to illustrate the Runge Phenomena
x = -1:0.5:1;
Yh = @(t) 1./(1+25*t.^2);
y = Yh(x);

F = genLagrangePolyInterp(x,y);

xMin = -1;
xMax = 1;
nX = 100;
xSpace = linspace(xMin,xMax,nX);
figure(2)
plot(xSpace,F(xSpace),...
    x,y,'sk',...
    xSpace,Yh(xSpace),'-r','linewidth',2,'markersize',10);
grid on
title('Lagrange Interpolation','fontsize',18,'fontweight','bold');
xlabel('X','fontsize',16,'fontweight','bold');
legend('Interpolated','Sample Points','True Function');
ylabel('Y','fontsize',16,'fontweight','bold');
set(gca,'fontsize',12,'fontweight','bold');

%% Compare with Chebychev Interpolation
%Add the Chebfun directories to path
Cheb_dir = '../../chebfun-master/chebfun-master';
addpath(Cheb_dir);

g = chebfun(Yh,[-1 1]);

figure(3)
plot(xSpace,g(xSpace),'-g',...
    xSpace,F(xSpace),...
    xSpace,Yh(xSpace),'--r','linewidth',2,'markersize',10);
grid on
title('Lagrange Interpolation','fontsize',18,'fontweight','bold');
xlabel('X','fontsize',16,'fontweight','bold');
legend('Chebychev Interp.','Lagrange Interp','True Function');
ylabel('Y','fontsize',16,'fontweight','bold');
set(gca,'fontsize',12,'fontweight','bold');
