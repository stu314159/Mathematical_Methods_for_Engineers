%% Lecture 33 MATLAB
clear
clc
close 'all'

%%
u_exact = @(x) x - sinh(x)/sinh(1);

%% Find a
Resid = @(x,a) -2*a - a.*x + a.*(x.^2) + x;
Weight = @(x) x - x.^2;
F = @(a) integral(@(x) Weight(x).*Resid(x,a),0,1);
a = fzero(F,1);
fprintf('a = %12.11f \n',a);
%%
u_galerkin = @(x) a*x.*(1-x);

a = 0; b = 1;
Nx = 100;
X = linspace(a,b,Nx);

figure(1)
plot(X,u_exact(X),'--r',...
    X,u_galerkin(X),'-b','linewidth',2)
title('Exact vs. Galerkin','fontsize',16,'fontweight','bold');
xlabel('X','fontsize',14,'fontweight','bold');
ylabel('u(X)','fontsize',14,'fontweight','bold');
grid on
set(gca,'fontsize',12,'fontweight','bold')
legend('Exact','Galerkin')