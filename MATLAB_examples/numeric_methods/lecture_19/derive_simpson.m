%% derive_simpson.m

%% does not yet work.  Problem is with use of fsolve and/or definition of simpsonF
%% Easy to do with EES.
% carry out non-linear solver routine necessary to derive Simpson's rule as
% at 3-point formula that will exactly integrate 3rd-order polynomials.
% This will be used to help explain why Simpson's Rule is a 3-point formula that
% achieves 4th order convergence as opposed to the Trapezoidal rule which
% is a 2-point formula that only achieves 2nd order convergence.

%%
% Also, this will be a good demonstration of the use of |fsolve| to handle
% a system of non-linear equations.

clear
clc
close all


%% Background
%
% Simpson's rule can be viewed as a 3-point integration scheme where we
% replace the function to be integrated -- f(x) -- with a 2nd order
% polynomial which we require to be equal to f(x) at the endpoints and the
% midpoint.  The resulting quadrature scheme can be said then to have the
% "sample points" x1 = a, x2 = (a+b)/2 and x3 = b.  The "weights" are 1/3,
% 4/3 and 1/3 respectively.  

% We derive the same rule if we relax the constraint on x2 and instead
% require that 3rd-order polynomials be computed exactly.  We use a
% starting estimate for x2 = 0.25 and set the weights initially to be all
% zero and assign them to the vector xW.  Feed this vector into fsolve
% along with the function |simpsonF.m| that applies the integration
% constraint (need a more full explanation on a handout.)

% naive initial guesses
xW = [0.25 .25 0.25 0.25];

w = fsolve(@simpsonF,xW);

fprintf('The middle sample point is: x = %g.\n',w(1));
fprintf('The weights are: \n');
format long
disp(w(2:end))
format short

