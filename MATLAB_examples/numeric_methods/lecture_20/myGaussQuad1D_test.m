% myGaussQuad1D_test.m

% clean the workspace
clear
clc
close all

% define function, interval and # Gauss Points
f = @(x) exp(-(x.^2));
a = -3;
b = 3;

P = 5;

% call myGaussQuad1D
intF = myGaussQuad1D(f,a,b,P);

% call built-in function integral
intF_builtIn = integral(f,a,b);


% compare the results.
fprintf('Result for my Gauss Quadrature for P = %i is: %18.16f \n',P,intF);
fprintf('Result for built in function integral = %18.16f \n',intF_builtIn);
fprintf('Relative difference = %g \n',abs(intF - intF_builtIn)/intF_builtIn);


%% Composite Gaussian Integration
% break the domain into N-1 subdomains.  Perform P-th order Gauss
% Quadrature on each subdomain and add up the results.  Do this for several
% values of P and plot the convergence results.

% this is a slow and loopy implementation, but it is easy to code and it
% works. In the future, this should be vectorized.

P_space = 2:2:10;
N_space = 1:10;

intArray = NaN(length(P_space),length(N_space));
relErrorArray = NaN(length(P_space),length(N_space));

% compute the integral over all values of N and all values of P
for n = 1:length(N_space)
    N = N_space(n);
    x_space = linspace(a,b,N+1);
    for p = 1:length(P_space)
        P = P_space(p);
        intArray(p,n) = myGaussQuad1Dcomposite(f,a,b,P,N);
        relErrorArray(p,n) = abs(intArray(p,n) - intF_builtIn)/intF_builtIn;
    end %p
end %n

% plot the results
figure
colorspec = {'b','r','g','c','m','y'};
for p = 1:length(P_space)
   loglog(N_space,relErrorArray(p,:),colorspec{p},'linewidth',2);
   if p==1
       hold on
   elseif p==length(P_space)
       hold off
   end
end
title('Convergence of Gauss Quadrature With Different P',...
    'fontsize',16,'fontweight','bold');
grid on
xlabel('Number of Intervals','fontsize',14,'fontweight','bold');
ylabel('Relative Error','fontsize',14,'fontweight','bold');
set(gca,'fontsize',12,'fontweight','bold');

