% mc_convergence.m

% This script is intended to show the rate of convergence for the
% Monte-Carlo type integration scheme.  Same problem as before, but the
% code is structured to allow a sequence of increasing number of sample
% points.  The relative error and sample size is then plotted on a log-log
% plot to show the convergence rate.
clear
clc
close all


N_min = 10;
N_max = 27; % <-- keep this less than 28 or so unless you have a lot of time on your hands.
N = N_min:N_max;

rel_error = NaN(length(N),1);

for i = N
   num_points = 2^i;
   fprintf('Estimating pi with %d sample points.\n',num_points);
   points = rand(num_points,2);
   frac = sum(sum(points .* points,2)<=1)/num_points;
   err_idx = i-N_min+1;
   rel_error(err_idx) = abs(pi - frac*4)/pi;
    
end

loglog(2.^N,rel_error,'linewidth',2)
title('Convergence rate for Monte Carlo estimation of Pi',...
    'fontweight','bold','fontsize',16);
xlabel('Number of Samples','fontweight','bold','fontsize',14);
ylabel('Relative Error','fontweight','bold','fontsize',14);
grid on

gage_line = 3.2*((2.^(N)).^(-0.5));
hold on
loglog(2.^N,gage_line,'-r','linewidth',2)
legend('Rel Err of MC Estimate','N^{-1/2} Convergence');

% notice that convergence is slow (~ N^(-1/2)) and not monotonic owing to
% the randomized nature of the algorithm.  Nonetheless for integrations
% with a very large number of degrees of freedom, MC methods are the
% methods of choice.