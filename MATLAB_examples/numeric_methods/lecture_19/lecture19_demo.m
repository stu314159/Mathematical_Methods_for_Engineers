%% Newton-Cotes Integration Demo
clear
clc
close 'all'

%% Setup
% define my function
f = @(x) exp(-x.^2);

% define the interval
xMin = -3;
xMax = 3;

n_exp = 5:9;
N = 2.^(n_exp) + 1; % since Simpson's needs an odd number of points.

% initialize data arrays
int_mp = nan(length(n_exp),1);
int_trap = nan(length(n_exp),1);
int_simp = nan(length(n_exp),1);
int_3p5o = nan(length(n_exp),1);
h_array = nan(length(n_exp),1);

err_mp = nan(length(n_exp),1);
err_trap = nan(length(n_exp),1);
err_simp = nan(length(n_exp),1);
err_3p5o = nan(length(n_exp),1);

% calculate "Gold Standard"
intF_GS = integral(f,xMin,xMax);

%% Do Calculations
% calculate integral with various methods
% and various numbers of subintervals.
h = nan(size(n_exp));
for s = 1:length(n_exp)
   h(s) = (xMax-xMin)/(N(s)-1);
   int_mp(s) = midpoint(f,xMin,xMax,N(s));
   err_mp(s) = norm(int_mp(s)-intF_GS,2)/norm(intF_GS,2);
   int_trap(s) = trapezoidal(f,xMin,xMax,N(s));
   err_trap(s) = norm(int_trap(s)-intF_GS,2)/norm(intF_GS,2);
   int_simp(s) = simpson(f,xMin,xMax,N(s));
   err_simp(s) = norm(int_simp(s)-intF_GS,2)/norm(intF_GS,2);
   int_3p5o(s) = ThreePointFifthOrder(f,xMin,xMax,N(s));
   err_3p5o(s) = norm(int_3p5o(s)-intF_GS,2)/norm(intF_GS,2);
end

%% Visualize Results
figure(1)
loglog(h,err_mp,'-r',...
    h,err_trap,'-c',...
    h,err_simp,'-b',...
    h,err_3p5o,'-g','linewidth',3)
grid on
xlabel('h','fontsize',12,'fontweight','bold');
ylabel('Relative Error','fontsize',12,'fontweight','bold');
title('Relative Error vs. h','fontsize',14,'fontweight','bold');
set(gca,'fontsize',10,'fontweight','bold');
legend('Midpoint','Trapezoidal','Simpson',...
    '3 point, 5th order','location','best');


%% Local Functions
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

function y =  ThreePointFifthOrder(f,xMin,xMax,N)
% in this context, N will be the number of subdivisions not
% the number of sample points.
xI = NaN(3,1);
xI(1) = 0.1127016653792583;
xI(2) = 0.5;
xI(3) = 0.887298334620755;

wI = NaN(3,1);
wI(1) = 0.277777777777777786;
wI(2) = 0.444444444444444493;
wI(3) = wI(1);

xSplit = linspace(xMin,xMax,N+1);
Jac = xSplit(2) - xSplit(1);

y = 0;
for d = 1:N
    a = xSplit(d); b = xSplit(d+1);
    xT = a + (b-a)*xI;
    y = y + f(xT)'*wI*Jac;
end

end
