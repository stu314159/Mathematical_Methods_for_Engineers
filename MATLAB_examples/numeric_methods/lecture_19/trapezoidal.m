function y = trapezoidal(f,xMin,xMax,N)
%trapezoidal(f,xMin,xMax,N) uses the trapezoidal rule with N equally spaced
%subintervals to estimate the integral of f over domain [xMin,xMax]
% inputs:
% f -- function_handle.  Handle to the function--f(x)-- to be integrated
% xMin -- scalar.  Lower bound of integration.
% xMax -- scalar. Upper bound of integration.
% N -- scalar.  Number of subdivisions.
% output:
% y -- scalar.  Approximate of the integral of f(x) from xMin to xMax.
%% create discrete spatial dimension with N subintervals
x = linspace(xMin,xMax,N+1); 
h = x(2)- x(1);

%% find f(x)
Fx = f(x); 

%% form the integration operator
intOp = h*ones(1,N+1); intOp(1)=0.5*h; intOp(end)=0.5*h;

%% apply the Trapezoidal integration operator
y = intOp*Fx';

end