% myGaussQuad1D_test.m
%% 
% clean the workspace
clear
clc
close all

% define function, interval and # Gauss Points
function_select = 1;

switch function_select
    
    case 1
        f = @(x) exp(-(x.^2));
        a = -3;
        b = 3;
        
    case 2
        f = @(x) (1./x).*exp(-x);
        a = 1;
        b = 100;
        
    case 3
        f = @(x) sinh(sqrt(2.29*x)).*exp(-x./0.965);
        a = 0;
        b = 100000000;
        
end



% define order of integration and number of subdivisions
P = 8;
N = 10000;

% call myGaussQuad1D
intF = gaussLegendre(f,a,b,P,N);

% call built-in function integral
intF_builtIn = integral(f,a,b);


% compare the results.
fprintf('Result for my Gauss Quadrature for P = %i is: %18.16f \n',P,intF);
fprintf('Result for built in function integral = %18.16f \n',intF_builtIn);
fprintf('Relative difference = %g \n',abs(intF - intF_builtIn)/intF_builtIn);


%% Convergence demonstration
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
        intArray(p,n) = gaussLegendre(f,a,b,P,N);
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