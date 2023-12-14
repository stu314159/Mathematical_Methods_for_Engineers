% legendre_expansion.m
% Purpose: expand a given function in terms of Legendre Polynomials and
% plot the results.
% Written by: CDR Stu Blair, USNA Mechanical Engineering Department
% Date: 27 July 2014

clear
clc
close all

% Set Maximum order of Legendre Polynomial for expansion
P = 8; 

% Pick a function
functionSelect = 5;

switch functionSelect
   
    case 1
        f = @(x) exp(-x.^2);
        a = -3; b = 3;
    case 2
        f = @(x) x.^2 - x.^3 + 3; % Only 0, 2 and 3 order terms non-zero
        a = -3; b = 3;
    case 3
        f = @(x) cos(x); % even function; odd order terms zero
        a = -4*pi; b = 4*pi;
    case 4
        f = @(x) sin(x); % odd function; even order terms zero
        a = -4*pi; b = 4*pi;
    case 5
        f = @(x) cosh(x) + sinh(x); % combination of both
        a = -2; b = 2;
    case 6
        f = @(x) square(x); % non-smooth function. *very* hard to approximate
        a = -4*pi; b = 4*pi;
    otherwise
        f = @(x) x; % default value
        a = -1; b = 1;
    
end

% generate Legendre Polynomials of order 0 through P
Pn = cell(P+1,1);
Pn{1} = @(x) 1;
Pn{2} = @(x) x;
if P == 0
    error('P must be greater than 0');
else
    for n = 2:P
        % use recurrence relation to generate higher order Legendre
        % Polynomials ("Pn functions")
        Pn{n+1} = @(x) (2*(n-1)+1)*x.*Pn{n}(x)./((n-1)+1) - (n-1)*Pn{n-1}(x)./((n-1)+1);
    end
end

% This is a change of variables from [-1,1] to [a,b]
% apply this to the Pn polynomials
xBack = @(t) (-2/(a-b))*t + (a+b)/(a-b);
JacBack= 2/(b-a);

% compute coefficients for Legendre Polynomial Exapnsion
Coef = NaN(P+1,1);
for n = 1:(P+1)
    fExp = @(x)f(x).*Pn{n}(xBack(x));
    Coef(n) = ((2*(n-1)+1)/2)*integral(fExp,a,b)*JacBack;
end
% Plot each term of the Legendre Polynomial Expansion as well as their sum.
colorspec = {'b','r','g','c','m','y'};
x = a:.01:b;

% Form Approximate:
errorNorm = NaN(P+1,1);
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,2,1)
fTot = ones(1,length(x))*Pn{1}(xBack(x))*Coef(1);
plot(x,fTot,'b','linewidth',2);
errorNorm(1) = norm((fTot - f(x)),2)/norm(f(x),2);
hold on
for p = 2:(P+1)
   fTerm = Coef(p)*Pn{p}(xBack(x));
   plot(x,fTerm,colorspec{mod(p,6)+1},'linewidth',2);
   fTot = fTot  + Coef(p)*Pn{p}(xBack(x));
   errorNorm(p) = norm((fTot - f(x)),2)/norm(f(x),2);
end
hold off
title('Legendre Polynomial Components of f(x)','fontweight','bold',...
    'fontsize',16);
xlabel('x','fontsize',14,'fontweight','bold');
ylabel('f(x)','fontsize',14,'fontweight','bold');

% make a plot of the coefficients
subplot(2,2,2)
bar(0:P,Coef)
%semilogy(0:P,abs(Coef),'linewidth',2);
title('Coefficients of the Legendre Polynomials','fontsize',16,...
    'fontweight','bold')
xlabel('Order of Legendre Polynomial Term','fontsize',14,'fontweight','bold');
ylabel('Value of Expansion Coefficient','fontsize',14,'fontweight','bold');

% plot the sum of all of the terms
subplot(2,2,3)
plot(x,fTot,'b','linewidth',2);
grid on
hold on
plot(x,f(x),'r','linewidth',2);
titleString = sprintf('Legendre approximate of f(x) with %i terms.',...
    P+1);
title(titleString,'fontweight','bold','fontsize',16);
xlabel('x','fontsize',14,'fontweight','bold');
ylabel('f(x)','fontsize',14,'fontweight','bold');
axis([min(x) max(x) min(f(x)) max(f(x))]);
legend('Approximate','Exact');

% plot the convergence behavior
subplot(2,2,4)
semilogy(0:P,errorNorm,'linewidth',2)
grid on
title('Convergence Behavior','fontweight','bold','fontsize',16);
xlabel('Order of Legendre Polynomial Approximation','fontsize',14,'fontweight','bold');
ylabel('Relative Error','fontsize',14,'fontweight','bold');
