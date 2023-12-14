%%midpoint.m
function y = midpoint(f,xMin,xMax,N)
% function y = midpoint(f,a,b,N)
% inputs:
% f -- function_handle.  Handle to the function--f(x)--to be integrated
% xMin -- scalar.  Lower bound of integration
% xMax -- scalar.  Upper bound of integration
% N -- scalar.  Number of subdivisions
% output:
% y -- scalar.  Approximate of the integral of f(x) from a to b.
xS = linspace(xMin,xMax,N+1);
xMid = (1/2)*(xS(1:(end-1))+xS(2:end));
h = xS(2)-xS(1);
y = h*sum(f(xMid));
end