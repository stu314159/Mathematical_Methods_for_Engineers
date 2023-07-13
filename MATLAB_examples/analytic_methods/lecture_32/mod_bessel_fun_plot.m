clear
clc
close 'all'

xMin = 1e-3;
xMax = 3;
Nx = 100;

X = linspace(xMin,xMax,Nx);

nu = 0;

figure(1)
plot(X,besseli(nu,X),'-b',...
    X,besselk(nu,X),'-g',...
    'LineWidth',3);
grid on
title('Modified Bessel Functions of Order Zero',...
    'FontSize',14,'FontWeight','bold');
xlabel('\alpha r','FontSize',12,'FontWeight','bold');
legend('I_0(\alpha r)','K_0(\alpha r)');
