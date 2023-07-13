%% script to plot J_1(x) and Y_1(x)

clear
clc
close 'all'

xMin = 1e-1; xMax = 15; Nx = 1000;
nu = 1;
X = linspace(xMin,xMax,Nx);
figure(1)
plot(X,besselj(nu,X),'-b',...
    X,bessely(nu,X),'-g',...
    'linewidth',3);
grid on
title('\textbf{Bessel functions of order }$$\nu=1$$',...
    'Interpreter','latex',...
    'FontSize',18,...
    'FontWeight','bold');
set(gca,'fontsize',12,'fontweight','bold');
xlabel('X','fontsize',14,'fontweight','bold');
legend('J_{1}(x)','Y_{1}(x)','Location','best');