% mc_convergence.m

% This script is intended to show the rate of convergence for the
% Monte-Carlo type integration scheme.  Same problem as before, but the
% code is structured to allow a sequence of increasing number of sample
% points.  The relative error and sample size is then plotted on a log-log
% plot to show the convergence rate.
clear
clc
close all
%% Parameters
f = @(x) exp(-x.^2);
a = -3; b = 3;
fMax = 1; fMin = 0;
% Compute area of bounding box
Abox = (b-a)*(fMax-fMin);

intGold = integral(f,a,b,'RelTol',1e-10);

%%

N_min = 10;
N_max = 28; % <-- keep this less than 28 or so unless you have a lot of time on your hands.
N = N_min:N_max;

rel_error = NaN(length(N),1);

for i = N
   num_points = 2^i;
   fprintf('Estimating integral with %d sample points.\n',num_points);
   x_s = a + (b-a)*rand(num_points,1);
   y_s = fMin + (fMax-fMin)*rand(num_points,1);
   frac = sum(y_s <= f(x_s))/num_points;
   err_idx = i-N_min+1;
   rel_error(err_idx) = abs(frac*Abox - intGold)/intGold;
    
end

%% Plot the results
loglog(2.^N,rel_error,'-b','linewidth',2)
title('Convergence rate for Hit or Miss',...
    'fontweight','bold','fontsize',16);
xlabel('Number of Samples','fontweight','bold','fontsize',14);
ylabel('Relative Error','fontweight','bold','fontsize',14);
grid on

gage_line = 3.2*((2.^(N)).^(-0.5));
hold on
loglog(2.^N,gage_line,'-r','linewidth',2)
legend('Rel Err of MC Estimate','N^{-1/2} Convergence');
axis([10^2 10^9 1e-5 0.25])
% notice that convergence is slow (~ N^(-1/2)) and not monotonic owing to
% the randomized nature of the algorithm.  Nonetheless for integrations
% with a very large number of degrees of freedom, MC methods are the
% methods of choice.