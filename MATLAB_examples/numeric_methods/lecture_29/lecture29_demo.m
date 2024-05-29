%% Lecture 29 Demo
clear
clc
close 'all'
%%
% 
%   Problem to be solved:
%
% $$\frac{d^2y}{dx^2}=y+sin(x) \ \ 0<x<2, \ \ y(0)=1, \ y(2) = 0$$
% 

C1 = (0.5*sin(2)-cosh(2))/sinh(2);
C2 = 1;
y_exact = @(x) C1*sinh(x) + C2*cosh(x) - 0.5*sin(x);

figure(1)
fplot(y_exact,[0,2],'linewidth',3);
title('Analytic Solution','fontsize',14,'fontweight','bold');
xlabel('X','fontsize',12,'fontweight','bold');
ylabel('Y(X)','fontsize',12,'fontweight','bold');
grid on

%% Solve with bvp4c
Ya = 1; Yb = 0; % specified boundary conditions
xMin = 0; xMax = 2; Nx = 5;
x = linspace(xMin,xMax,Nx);
Yguess = [1 0];
F = @(x,w) [w(2); w(1)+sin(x)]; 
bcfun = @(ya,yb) [ya(1)- Ya; yb(1)-Yb];
%solinit = bvpinit(x,Yguess);
solinit = bvpinit([xMin xMax],Yguess);
options = bvpset('RelTol',1e-10,'AbsTol',1e-8,'NMax',5000);
sol = bvp4c(F,bcfun,solinit,options);

figure(2)
plot(sol.x,sol.y(1,:),'-sb','linewidth',3);
title('Solution with BVP4C','fontsize',14,...
    'fontweight','bold');
xlabel('X','fontsize',12,'fontweight','bold');
ylabel('Y(X)','fontsize',12,'fontweight','bold');
grid on

% get error
rel_err = norm(y_exact(sol.x)-sol.y(1,:),2)/...
    norm(y_exact(sol.x),2);
fprintf('Relative error = %g \n',rel_err);
fprintf('Number of grid-points used: %u \n',...
    length(sol.x));

%% Solve with BVP5C
sol5c = bvp5c(F,bcfun,solinit,options);
figure(3)
plot(sol5c.x,sol5c.y(1,:),'-sb','linewidth',3);
title('Solution with BVP5C','fontsize',14,...
    'fontweight','bold');
xlabel('X','fontsize',12,'fontweight','bold');
ylabel('Y(X)','fontsize',12,'fontweight','bold');
grid on

% get error
rel_err = norm(y_exact(sol5c.x)-sol5c.y(1,:),2)/...
    norm(y_exact(sol5c.x),2);
fprintf('Relative error = %g \n',rel_err);
fprintf('Number of grid-points used: %u \n',...
    length(sol5c.x));

%% Example 2
F = @(x,t) ex2(x,t);
L = 4e-3; % m, length of wire
xMin = 0; xMax = L/2; N = 500;
x = linspace(xMin,xMax,N);

To = 300; % K, temperature at x=0
bcfun = @(ya,yb) [ya(1)- To; yb(2) - 0];

Tguess = [To 0]; % second entry is dT/dx
solinit = bvpinit(x,Tguess);
options = bvpset('RelTol',1e-5,'AbsTol',1e-5,'NMax',5000);
sol2 = bvp5c(F,bcfun,solinit,options);

figure(4)
plot(sol2.x,sol2.y(1,:),'-b','linewidth',3);
grid on
title('Example 2 Solution','fontsize',14,...
    'fontweight','bold');
xlabel('X [m]','fontsize',12,'fontweight','bold');
ylabel('T(X) [K]','fontsize',12,'fontweight','bold');
set(gca,'fontsize',10,'fontweight','bold');


%% Local Functions
function dTdx = ex2(~,T)
k = 72; % W/m-K
h = 2000; % W/m^2-K
emiss = 0.1;
sigma_sb = 5.67e-8; % units...
I = 2; % amps
rho_e = 32e-8; % Ohm-m
Tinf = 300; % K
D = 7.62e-5; % m

% make some terms to simplify the equation
c1 = 4*h/(k*D);
c2 = 4*emiss*sigma_sb/(k*D);
c3 = -I^2*rho_e/(k*(0.25*pi*D^2)^2);
dTdx = [T(2);
    c1*(T(1) - Tinf) + c2*(T(1).^4 - Tinf^4) + c3];

end