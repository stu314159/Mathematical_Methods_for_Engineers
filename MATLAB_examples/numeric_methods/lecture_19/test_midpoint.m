%% Numeric Integration: Midpoint Rule
% 

clear
clc
close all

% define my function
f = @(x) exp(-x.^2);


% define the interval
xMin = -3;
xMax = 3;

% Specify number of sub-divisions
N = 100;

% call the solver
intF = midpoint(f,xMin,xMax,N);

%intF_Analytic = 1/3; % for this case, I know the analytic solution.
intF_Analytic = integral(f,xMin,xMax);

% compute relative error

relError = abs(intF - intF_Analytic)/intF_Analytic;

fprintf('Integration complete! intF = %g.  Relative Error = %g.\n',...
    intF,relError);


%% Convergence Rate of the Midpoint Method

N = 1000;
intF = midpoint(f,xMin,xMax,N);
relError = abs(intF - intF_Analytic)/intF_Analytic;
fprintf('Integration complete! intF = %g.  Relative Error = %g.\n',...
    intF,relError);

%%
% As you can see, we increased the number of sub-divisions by a factor of
% 10, and the accuracy increased by a factor of 100.  Not a bad bargain.  A
% more detailed analysis shows that the relative error decreases as a
% function of the interval size squared.  In mathematical terms, this is
% referred to as "second order convergence".
Nspace = 2.^(7:15);
intF_vec = NaN(1,length(Nspace));
relError_vec = intF_vec;
for n = 1:length(Nspace)
    intF_vec(n) = midpoint(f,xMin,xMax,Nspace(n));
    relError_vec(n) = abs(intF_vec(n) - intF_Analytic)/intF_Analytic;
end
loglog(Nspace,relError_vec,'-b','linewidth',2);
grid on
hold on
loglog(Nspace,.02*(1./Nspace).^2,'--r','linewidth',2)
title('Convergence of Midpoint Rule','fontsize',16,'fontweight','bold');
xlabel('Number of Subdivisions','fontsize',14,'fontweight','bold');
ylabel('Relative Error','fontsize',14,'fontweight','bold');
legend('Relative Error','Second Order Convergence');
set(gca,'fontweight','bold');
%% Exact Integration of Linear Functions
% 
f = @(x) (pi/4).*x;
N = 10;
intF = midpoint(f,xMin,xMax,N);
intF_Analytic = (pi/8);
relError = abs(intF - intF_Analytic)/intF_Analytic;
fprintf('Integration complete! intF = %g.  Relative Error = %g.\n',...
    intF,relError);
%%
% 




