%% Assignment #6 Solution
clear
clc
close 'all'

%% Problem #1 (based on 6-34)
T = [50 100 150 200 400 600 800 1000]'; % K, silicon temperature
k = [28 9.1 4.0 2.7 1.1 0.6 0.4 0.3]'; % W/m-K, thermal conductivity

% jut plot the data
figure(1)
plot(T,k,'k*')
title('Thermal Conductivity vs. Temperature',...
    'fontsize',14,'fontweight','bold');
xlabel('Temperature [K]','fontsize',12,'fontweight','bold');
ylabel('Thermal Conductivity [W/m-K]','fontsize',12,...
    'fontweight','bold');
grid on
%%
% a)y = bx^m;  ln(y) = ln(b) + m*ln(x)

Y = log(k);
X = [T.^0 log(T)];
p = X'*Y;
c = (X'*X)\p;
% c(1) = ln(b) --> b = exp(c(1))
% c(2) = m
b = exp(c(1));
m = c(2);

fita = @(T) b*T.^m;

figure(2)
plot(T,k,'k*',...
    T,fita(T),'-c','linewidth',3,'markersize',10)
title('Thermal Conductivity vs. Temperature',...
    'fontsize',14,'fontweight','bold');
xlabel('Temperature [K]','fontsize',12,'fontweight','bold');
ylabel('Thermal Conductivity [W/m-K]','fontsize',12,...
    'fontweight','bold');
grid on

resid_a = k - fita(T);
err_a = resid_a'*resid_a;
fprintf('error for fit a: %g \n',err_a);
%%
% b) y = b*exp(m*x) --> log(y) = log(b) + m*x

X = [T.^0 T];
c = (X'*X)\(X'*Y); %<- same Y as before

% c(1) = log(b) --> b = exp(c(1))
% c(2) = m
b = exp(c(1)); m = c(2);
fitb = @(T) b*exp(m*T);

figure(3)
plot(T,k,'k*',...
    T,fitb(T),'-c')
title('Thermal Conductivity vs. Temperature',...
    'fontsize',14,'fontweight','bold');
xlabel('Temperature [K]','fontsize',12,'fontweight','bold');
ylabel('Thermal Conductivity [W/m-K]','fontsize',12,...
    'fontweight','bold');
grid on

resid_b = k - fitb(T);
err_b = resid_b'*resid_b;
fprintf('error for fit b: %g \n',err_b);
%%
% c) y = b*10^(mx) --> log10(y) = log10(b) + m*x
Y = log10(k);
X = [T.^0 T.^1];
p = X'*Y;
c = (X'*X)\p;
% c(1) = log10(b) --> b = 10^c(1), m = c(2)
b = 10^(c(1)); m = c(2);
fitc = @(T) b*10.^(m*T);

figure(4)
plot(T,k,'k*',...
    T,fitc(T),'-c')
title('Thermal Conductivity vs. Temperature',...
    'fontsize',14,'fontweight','bold');
xlabel('Temperature [K]','fontsize',12,'fontweight','bold');
ylabel('Thermal Conductivity [W/m-K]','fontsize',12,...
    'fontweight','bold');
grid on

resid_c = k - fitc(T);
err_c = resid_c'*resid_c;
fprintf('error for fit c: %g \n',err_c);
%%
% d) y = 1/(mx+b) --> 1/y = mx + b
Y = 1./k;
X = [T.^0 T.^1];
p = X'*Y;
c = (X'*X)\p;
b = c(1); m = c(2);
fitd = @(T) 1./(m*T + b);

T2 = linspace(min(T),max(T),1000);
figure(5)
plot(T,k,'k*',...
    T,fitd(T),'-c')
title('Thermal Conductivity vs. Temperature',...
    'fontsize',14,'fontweight','bold');
xlabel('Temperature [K]','fontsize',12,'fontweight','bold');
ylabel('Thermal Conductivity [W/m-K]','fontsize',12,...
    'fontweight','bold');
grid on

resid_d = k - fitd(T);
err_d = resid_d'*resid_d;
fprintf('error for fit d: %g \n',err_d);
%%
% e)y = mx/(b + x) --> Y = 1/y, X = 1/x, c1 = b/m, c2 = 1/m
% --> m = 1/c2; b = c1*m
Y = 1./k;
X = [T.^-1 T.^0];
p = X'*Y;
c = (X'*X)\p;
m = 1./c(2);
b = c(1)*m;

fite = @(x) m*x./(b + x);
figure(6)
plot(T,k,'k*',...
    T,fite(T),'-c')
title('Thermal Conductivity vs. Temperature',...
    'fontsize',14,'fontweight','bold');
xlabel('Temperature [K]','fontsize',12,'fontweight','bold');
ylabel('Thermal Conductivity [W/m-K]','fontsize',12,...
    'fontweight','bold');
grid on

resid_e = k - fite(T);
err_e = resid_e'*resid_e;
fprintf('error for fit e: %g \n',err_e);
%% Problem #3
T = (1e3)*[5 7.5 10 12.5 15 17.5 20 22.5 25 27.5 30]';
h = [3.3 7.5 41.8 51.8 61 101.1 132.9 145.5 171.4 225.8 260.9]';

T_int = 5000:500:3e4;
F_lagrange = genLagrangePolyInterp(T,h);

figure(7)
plot(T,h,'k*',...
    T_int,F_lagrange(T_int),'-c','linewidth',3,'markersize',8);
grid on
title('Problem #3 - Lagrange Polynomial Interpolation',...
    'fontsize',14,'fontweight','bold');
xlabel('T [K]','fontsize',12,'fontweight','bold');
ylabel('h [MJ/kg]','fontsize',12,'fontweight','bold');

method = 'spline';
F_spline = @(xi) interp1(T,h,xi,method);
figure(8)
plot(T,h,'k*',...
    T_int,F_spline(T_int),'-c','linewidth',3,'markersize',8);
grid on
title('Problem #3 - MATLAB Spline Interpolation',...
    'fontsize',14,'fontweight','bold');
xlabel('T [K]','fontsize',12,'fontweight','bold');
ylabel('h [MJ/kg]','fontsize',12,'fontweight','bold');

%% Local function for Lagrange Polynomial
function F = genLagrangePolyInterp(X,Y)
% function F = genLagrangePoly(X,Y) generates a Lagrange Polynomial that
% may be used to interpolate a function
% Inputs
% X = x-values of a function
% Y = f(X) for some function
%
% Outputs
% F - a function handle with the Lagrange interpolant

n = length(X); 

F = @(x) 0; %<-- start with zero

for i = 1:n
    L = @(x) Y(i); %<-- initialize the Lagrange Function
    for j = 1:n
        if j ~= i
            L = @(x) L(x).*((x - X(j))./(X(i) - X(j)));
        end
    end
    F = @(x) F(x)+L(x);
end
end