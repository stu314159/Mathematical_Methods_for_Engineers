%% Numeric Integration: Simpson's Rule
% In many respects, Simpson's Rule is similar to the Trapezoidal Rule.
% Sometimes call the "Simpson's 1/3 Rule", this method requires that there
% be an *even* number of equally sized intervals.  
%%
% The approximate value of the integral is given as:
%%
%
% $$ I(f) \approx \frac{h}{3} \left[ f(a) + 4\sum_{i=2,4,6,...}^{N} + 2\sum_{j = 3,5,7,...}^{N-1} + f(b) \right] $$
%
%
%%
% The function |simpson.m| is implemented as follows:
%%
%
%   function y = simpson(f,xMin,xMax,N)
% 
%   % create discrete spatial dimension with N subintervals
%   % make sure that N is even
%   if(mod(N,2)~=0)
%       N = N+1; % make the number of subintervals even
%   end
% 
%   x = linspace(xMin,xMax,N+1);
%   h = x(2) - x(1);
% 
%   % find f(x)
%   Fx = f(x);
% 
%   % form the integration operator
%   intOp = (h/3)*ones(1,N+1);
%   intOp(2:2:end) = 4*intOp(2:2:end);
%   intOp(3:2:end-1) = 2*intOp(3:2:end-1);
% 
%   % apply the Simpson's integration operator
% 
%   y = intOp*Fx';
% 
%   end
%
%%
% Notice the similarities between the implementation of the Simpson's
% Method and the Trapezoidal Method.
%%
% We will test the method in a similar way as well:
%%
clear
clc
close all

% define my function
%f = @(x) x.^4;
f = @(x) exp(-x.^2);

% define the interval
% xMin = 0;
% xMax = 1;
xMin = -3;
xMax = 3;

% Specify number of sub-divisions
N = 100;

% call the solver
intF = simpson(f,xMin,xMax,N);

%intF_Analytic = 1/5; % for this case, I know the analytic solution.
intF_Analytic = integral(f,xMin,xMax);

% compute relative error

relError = abs(intF - intF_Analytic)/intF_Analytic;

fprintf('Integration complete! intF = %g.  Relative Error = %g.\n',...
    intF,relError);

%% Convergence Rate of Simpson's Rule
% As with the Trapezoidal Method, with Simpson's Rule, if you would like a
% more accurate estimate of the definite integral, all you need to do is
% specify a larger number of sub-intervals. (and have your *computer* do
% correspondingly more work!)
%
Nspace = 2.^(5:12);
intF_vec = NaN(1,length(Nspace));
relError_vec = intF_vec;
for n = 1:length(Nspace)
    intF_vec(n) = simpson(f,xMin,xMax,Nspace(n));
    relError_vec(n) = abs(intF_vec(n) - intF_Analytic)/intF_Analytic;
end
loglog(Nspace,relError_vec,'-b','linewidth',2);
grid on
hold on
loglog(Nspace,.02*(1./Nspace).^4,'--r','linewidth',2)
title('Convergence of Simpsons Rule','fontsize',16,'fontweight','bold');
xlabel('Number of Subdivisions','fontsize',14,'fontweight','bold');
ylabel('Relative Error','fontsize',14,'fontweight','bold');
legend('Relative Error','Fourth Order Convergence');
set(gca,'fontweight','bold');
%%
% Notice that we achieve fourth-order convergence with Simpson's Rule: if
% we reduce our interval size by a factor of 2, our relative error will be
% reduced by a factor of 2^4 = 16.  