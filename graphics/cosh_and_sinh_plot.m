% cosh and sinh plot
clear
clc
close 'all'

xMin = -2;
xMax = 3;
Nx = 1000;
X = linspace(xMin,xMax,Nx);

figure (1)
plot(X,cosh(X),'-c^',...
    X,sinh(X),'-ks',...
    'linewidth',3,...
    'MarkerIndices',1:100:Nx)
xline(0,'linewidth',4);
yline(0,'linewidth',4);
grid on
%title('COSH and SINH','fontsize',14,'fontweight','bold');
xlabel('X')
legend('cosh(x)','sinh(x)','fontsize',12)