%% Lecture 30 MATLAB Demo
clear
clc
close 'all'

%% Example 1
R = 1.5e-2; % m, radius of fuel
w = 3.0e-3; % m, thickness of cladding
k = 16.75; % W/(m-K), thermal conductivity of clad
Q = 1e8; % W/m^2, Source term from heat dep in cladding.

Q2 = 6.32e5; % W/m^2, Heat flux due to heat produced in fuel.
T_inf = 423; % K, temperature of water flowing on cladding
h = 1e4; % W/(m^2-K), convective heat transfer coefficient

F = @(r,T) [T(2); -(1./r)*T(2) - Q./(k*r)*exp(-r/R)];
bcfun = @(Ta,Tb) [Ta(2)+ Q2/k; ...
    Tb(2) + (h/k)*(Tb(1)-T_inf)];
rMin = R; rMax = R+w; nR = 200;
Tguess = [T_inf 0];
solinit = bvpinit(linspace(rMin,rMax,nR),Tguess);

sol = bvp5c(F,bcfun,solinit);

figure(1)
plot(sol.x,sol.y(1,:),'-b','linewidth',2);
title('Fuel Clad Temperature','fontsize',14,...
    'fontweight','bold');
xlabel('R [m]','fontsize',12,'fontweight','bold');
ylabel('T(R) [K]','fontsize',12,'fontweight','bold');
grid on
set(gca,'fontsize',10,'fontweight','bold');

%% Parameter Example
%%
% 
%  Suppose that Q could be modified.  We migth conceive of a problem where
%  we want to find the value of Q such that the T(r=R) hits a specified
%  value.  MATLAB's built-in solvers provide functionality that allows this
%  type of analysis.
% 
R = 1.5e-2; % m, radius of fuel
w = 3.0e-3; % m, thickness of cladding
k = 16.75; % W/(m-K), thermal conductivity of clad
%Q = 1e8; % W/m^2, Source term from heat dep in cladding. 
Q2 = 6.32e5; % W/m^2, Heat flux due to heat produced in fuel.
T_inf = 423; % K, temperature of water flowing on cladding
h = 1e4; % W/(m^2-K), convective heat transfer coefficient

TgtTemp = 1500; % K

% add Q to my arguments for F
F = @(r,T,Q) [T(2); -(1./r)*T(2) - Q./(k*r)*exp(-r/R)];

% ... and for my BC functions
bcfun = @(Ta,Tb,Q) [Ta(2)+ Q2/k; ...
    Tb(2) + (h/k)*(Tb(1)-T_inf); ...
    Ta(1) - TgtTemp];
rMin = R; rMax = R+w; nR = 200;
Tguess = [T_inf 0];

% ... and add this to solinit
Qguess = 1e8;
solinit = bvpinit(linspace(rMin,rMax,nR),Tguess,Qguess);

sol2 = bvp5c(F,bcfun,solinit);

figure(2)
plot(sol2.x,sol2.y(1,:),'-b','linewidth',2);
title('Fuel Clad Temperature','fontsize',14,...
    'fontweight','bold');
xlabel('R [m]','fontsize',12,'fontweight','bold');
ylabel('T(R) [K]','fontsize',12,'fontweight','bold');
grid on
set(gca,'fontsize',10,'fontweight','bold');

fprintf('Q = %g \n',sol2.parameters);
