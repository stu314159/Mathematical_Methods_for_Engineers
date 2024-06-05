%% Lecture 33 MATLAB
clear
clc
close 'all'

%%
u_exact = @(x) x - sinh(x)/sinh(1);

u_trial = @(x,a) a.*(x - x.^2);
%% Compute the Residual

Resid = @(x,a) -2*a - a.*x + a.*(x.^2) + x;

%% Method of Least Squares

Weight_LS = @(x) -2 - (x - x.^2);
F_LS = @(a) integral(@(x) Weight_LS(x).*Resid(x,a),0,1);
a_LS = fzero(F_LS,1);

fprintf('a_LS = %12.11f \n',a_LS);

%% Co-location method
x_i = 0.5; % colocation point
F_CL = @(a) Resid(x_i,a);

a_CL = fzero(F_CL,1);

fprintf('a_CL = %12.11f \n',a_CL);

%% Galerkin Method

Weight_Galerkin = @(x) x - x.^2;
F = @(a) integral(@(x) Weight_Galerkin(x).*Resid(x,a),0,1);
a = fzero(F,1);
fprintf('a = %12.11f \n',a);
%% Plot Solution
u_CL = @(x) u_trial(x,a_CL);
u_LS = @(x) u_trial(x,a_LS);
u_galerkin = @(x) u_trial(x,a);

a = 0; b = 1;
Nx = 100;
X = linspace(a,b,Nx);

figure(1)
plot(X,u_exact(X),'--r',...
    X, u_LS(X),'-k',...
    X, u_CL(X),'-g',...)
    X,u_galerkin(X),'-b','linewidth',2)
title('Exact vs. MWR','fontsize',16,'fontweight','bold');
xlabel('X','fontsize',14,'fontweight','bold');
ylabel('u(X)','fontsize',14,'fontweight','bold');
grid on
set(gca,'fontsize',12,'fontweight','bold')
legend('Exact','Least Squares','Co-location','Galerkin')

%% Weak Form
W = @(x) x - x.^2;
Ut = @(x,a) a.*(x - x.^2);
dW_dx = @(x) 1 - 2*x;
dUt_dx = @(x,a) a.*(1-2*x);

Weak_Form = @(a) W(1)*dUt_dx(1,a) - W(0)*dUt_dx(0,a) - ... % boundary term
    integral(@(x) dUt_dx(x,a).*dW_dx(x),0,1) - ...
    integral(@(x) W(x).*Ut(x,a),0,1) + ...
    integral(@(x) W(x).*x,0,1);

a_weak = fzero(Weak_Form,1);

fprintf("Weak Form solution, a = %12.11f \n",a_weak);