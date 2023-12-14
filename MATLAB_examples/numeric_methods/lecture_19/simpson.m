function y = simpson(f,xMin,xMax,N)
%simpson(f,xMin,xMax,N) uses Simpson's rule to approximate the integral of
%f over domain [xMin,xMax].
% inputs:
% f -- function_handle.  Handle to the function--f(x)-- to be integrated
% xMin -- scalar.  Lower bound of integration.
% xMax -- scalar. Upper bound of integration.
% N -- scalar.  Number of subdivisions.
% output:
% y -- scalar.  Approximate of the integral of f(x) from xMin to xMax.

%% create discrete spatial dimension with N subintervals
% make sure that N is even
if(mod(N,2)~=0)
    N = N+1; % make the number of subintervals even
end

x = linspace(xMin,xMax,N+1);
h = x(2) - x(1);

%% find f(x)
Fx = f(x);

%% form the integration operator
intOp = (h/3)*ones(1,N+1);
intOp(2:2:end) = 4*intOp(2:2:end);
intOp(3:2:end-1) = 2*intOp(3:2:end-1);

%% apply the Simpson's integration operator
y = intOp*Fx';

end