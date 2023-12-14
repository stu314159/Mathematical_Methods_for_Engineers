%% Numeric Integration: Trapezoidal Rule
% Probably the simplest method for numeric integration is the Trapezoidal
% Rule.  The form of Trapezoidal Rule illustrated here is called the
% Composite Trapezoidal method.  It requires that the domain of integration
% be discretized into subintervals.  The function to be
% integrated is evaluated at the end-points of each sub-interval and
% employs a linear approximation for the function between the end-points.
% Note that there is no requirement that each subinterval be the same size.
%  Nonetheless, we will also, for this analysis, use equally spaced
%  sub-intervals; it is slightly easier to implement that way and its
%  convergence properties are far more easy to analyze.

%%
% 
% <<../TrapRule.gif>>
% 

%%
% The approximate value of the integral is:
%%
% $$ I(f) \approx \frac{h}{2} \left[ f(a) + f(b) \right] + h\sum_{i =
% 2}^{N} f\left(x_{i}\right)$$
%%
%
%
%%
% The function |trapezoidal.m| is implemented as follows:
%%
%
%   function y = trapezoidal(f,xMin,xMax,N)
% 
%   % create discrete spatial dimension with N subintervals
%   x = linspace(xMin,xMax,N+1); 
%   h = x(2)- x(1);
% 
%   % find f(x)
%   Fx = f(x); 
% 
%   % form the integration operator
%   intOp = h*ones(1,N+1); intOp(1)=0.5*h; intOp(end)=0.5*h;
% 
%   % apply the Trapezoidal integration operator
%   y = intOp*Fx';
% 
%   end
%
%%
% To test this method, we will define a function and interval over which it
% will be integrated.  To be sure that the method is implemented correctly,
% we will chose a function for which we can find the definite integral
% exactly.  For this example, we choose:
%%
% $$ f(x) = x^{2} \ \ x \ \in \ [0,1] $$
% 
%%
% We can easily calculate the integral of this function; it is equal to
% 1/3.

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
intF = trapezoidal(f,xMin,xMax,N);

intF_Analytic = integral(f,xMin,xMax);

% compute relative error

relError = abs(intF - intF_Analytic)/intF_Analytic;

fprintf('Integration complete! intF = %g.  Relative Error = %g.\n',...
    intF,relError);


%% Convergence Rate of the Trapezoidal Method
% It is probably intuitively clear to you that if you want a more exact
% value of the integral, you merely need to increas the number of
% sub-intervals.
N = 1000;
intF = trapezoidal(f,xMin,xMax,N);
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
    intF_vec(n) = trapezoidal(f,xMin,xMax,Nspace(n));
    relError_vec(n) = abs(intF_vec(n) - intF_Analytic)/intF_Analytic;
end
loglog(Nspace,relError_vec,'-b','linewidth',2);
grid on
hold on
loglog(Nspace,.02*(1./Nspace).^2,'--r','linewidth',2)
title('Convergence of Trapezoidal Rule','fontsize',16,'fontweight','bold');
xlabel('Number of Subdivisions','fontsize',14,'fontweight','bold');
ylabel('Relative Error','fontsize',14,'fontweight','bold');
legend('Relative Error','Second Order Convergence');
set(gca,'fontweight','bold');
%% Exact Integration of Linear Functions
% Since the Trapezoidal Rule approximates the function in a piece-wise
% linear way, any linear function will be integrated *exactly*.
f = @(x) (pi/4).*x;
N = 10;
intF = trapezoidal(f,xMin,xMax,N);
intF_Analytic = (pi/8);
relError = abs(intF - intF_Analytic)/intF_Analytic;
fprintf('Integration complete! intF = %g.  Relative Error = %g.\n',...
    intF,relError);
%%
% Notice that though I indicate that the integration is exact, there is
% still a non-zero relative error.  Why is this? _The result is exact up to
% the relative precision of double precision floating point arithmetic.
% Effectively, since the relative error is less than |eps|, from the
% perspective of floating point number representatation, there is no
% difference between the given relative error and 0_.




