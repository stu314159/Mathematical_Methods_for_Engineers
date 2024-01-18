%% RK-Method Convergence Tester
clear 
clc
close 'all'

%% Select Solver
Solver_choice = 13;

% Explicit Methods
% 1 = explicit Euler
% 2 = 2nd order RK - Heun's Method (10.5.1)
% 3 = 3rd order RK - "Classical third-order RK" (10.5.2)
% 4 = 3rd order RK - "Nystrom Scheme" (Iserles 3.2 pg 40-41)
% 5 = 2/3 order RK - "Bogacki-Shampine method" (see Wiki - used for ode23
% for MATLAB)
% 6 = 4th order RK - "Classical fourth-order RK" (10.5.3)
% 7 = 4/5th order RK - "Fehlberg Scheme" (see Wikipedia - default ode45 for Octave)
% 8 = 4/5th order RK - "Dormand-Prince Scheme" (see Wiki - default ode45 for MATLAB)

% Implicit Methods
% 9 = 2-stage, 3rd order implicit RK
% 10 = 2-stage, 4th order implicit RK
% 11 = 4th order Lobatto IIIC (try for stiff problems)
% 12 = 1 stage Backward Euler
% 13 = implicit midpoint
% 14 = Gauss Legendre method of order 6
% 15 = 4th order IRK from Butcher, Eq 237c
% 16 - Lobatto III (5-stage, 8th-order), Butcher pg 228

switch Solver_choice
    
    case 1
        solver = @odesEULER;
        BT = [];
        order = 1;
        
    case 2
        solver = @odesExplicitRK;
        s = 2;
        BT = zeros(s+1,s+1);
        C = [0; 2/3; 0];
        B = [0 1/4 3/4];
        A = [0 0;
            2/3 0;];
        BT(:,1) = C;
        BT(end,:) = B;
        BT(1:s,2:end) = A;
        order = 2;
        
    case 3
        solver = @odesExplicitRK;
        s = 3;
        BT = zeros(s+1,s+1);
        C = [0; 1/2; 1; 0];
        B = [0 1/6 4/6 1/6];
        A = [0 0 0;
            1/2 0 0;
            -1 2 0;];
        
        BT(:,1) = C;
        BT(end,:) = B;
        BT(1:s,2:end) = A;
        order = 3;
        
    case 4
        solver = @odesExplicitRK;
        s = 3;
        BT = zeros(s+1,s+1);
        C = [0; 2/3; 2/3; 0];
        B = [0 1/4 3/8 3/8];
        A = [0 0 0;
            2/3 0 0;
            0 2/3 0;];
        
        BT(:,1) = C;
        BT(end,:) = B;
        BT(1:s,2:end) = A;
        order = 3;
    case 5
        solver = @odesExplicitRK;
        s = 4;
        BT = zeros(s+1,s+1);
        C = [0; 1/2; 3/4; 1; 0];
        B = [0 2/9 1/3 4/9 0];
        A = [0 0 0 0;
            1/2 0 0 0;
            0 3/4 0 0;
            2/9 1/3 4/9 0];
        BT(:,1) = C;
        BT(end,:) = B;
        BT(1:s,2:end) = A;
        order = 3;
        
    case 6
        solver = @odesExplicitRK;
        s = 4;
        BT = zeros(s+1,s+1);
        C = [0; 1/2; 1/2; 1; 0];
        B = [0 1/6 2/6 2/6 1/6];
        A = [0 0 0 0;
            1/2 0 0 0;
            0 1/2 0 0;
            0 0 1 0;];
        BT(:,1) = C;
        BT(end,:) = B;
        BT(1:s,2:end) = A;
        order = 4;
        
    case 7
        solver = @odesExplicitRK;
        s = 6;
        BT = zeros(s+1,s+1);
        C = [0; 1/4; 3/8; 12/13; 1; 1/2; 0];
        B = [0 16/135 0 6656/12825 28561/56430 -9/50 2/55];
        A = [0 0 0 0 0 0;
            1/4 0 0 0 0 0;
            3/32 9/32 0 0 0 0;
            1932/2197 -7200/2197 7296/2197 0 0 0;
            439/216 -8 3680/513 -845/4104 0 0;
            -8/27 2 -3544/2565 1859/4104 -11/40 0;];
        BT(:,1) = C;
        BT(end,:) = B;
        BT(1:s,2:end) = A;
        order = 5;
        
    case 8
        solver = @odesExplicitRK;
        s = 7;
        BT = zeros(s+1,s+1);
        C = [0; 1/5; 3/10; 4/5; 8/9; 1; 1; 0];
        B = [0 35/384 0 500/1113 125/192 -2187/6784 11/84 0];
        A = [0 0 0 0 0 0 0;
            1/5 0 0 0 0 0 0;
            3/40 9/40 0 0 0 0 0;
            44/45 -56/15 32/9 0 0 0 0;
            19372/6561 -25360/2187 64448/6561 -212/729 0 0 0;
            9017/3168 -355/33 46732/5247 49/176 -5103/18656 0 0;
            35/384 0 500/1113 125/192 -2187/6784 11/84 0;];
        BT(:,1) = C;
        BT(end,:) = B;
        BT(1:s,2:end) = A;
        order = 5;
        
    case 9
        solver = @odesImplicitRK;
        s = 2;
        BT = zeros(s+1,s+1);
        C = [0; 2/3; 0];
        B = [0 1/4 3/4];
        A = [1/4 -1/4;
            1/4 5/12];
        BT(:,1) = C;
        BT(end,:) = B;
        BT(1:s,2:end) = A;
        order = 3;
        
    case 10
        solver = @odesImplicitRK;
        s = 2;
        BT = zeros(s+1,s+1);
        C = [1/2-sqrt(3)/6; 1/2+sqrt(3)/6; 0];
        B = [0, 1/2, 1/2];
        A = [1/4 1/4-sqrt(3)/6;
            1/4+sqrt(3)/6, 1/4];
        BT(:,1) = C;
        BT(end,:) = B;
        BT(1:s,2:end) = A;
        order = 4;
        
    case 11
        solver = @odesImplicitRK;
        s = 3;
        BT = zeros(s+1,s+1);
        C = [0; 0.5; 1; 0];
        B = [0, 1/6, 2/3, 1/6];
        A = [1/6, -1/3, 1/6;
            1/6, 5/12, -1/12;
            1/6, 2/3, 1/6];
        
        BT(:,1) = C;
        BT(end,:) = B;
        BT(1:s,2:end) = A;
        order = 4;
        
    case 12
        solver = @odesImplicitRK;
        s = 1;
        BT = zeros(s+1,s+1);
        C = [1; 0];
        B = [0, 1];
        A = 1;
        BT(:,1) = C;
        BT(end,:) = B;
        BT(1:s,2:end) = A;
        order = 1;
        
    case 13
        solver = @odesImplicitRK;
        s = 1;
        BT = zeros(s+1,s+1);
        C = [1/2; 0];
        B = [0,1];
        A = 1/2;
        BT(:,1) = C;
        BT(end,:) = B;
        BT(1:s,2:end) = A;
        order = 2;
        
    case 14
        solver = @odesImplicitRK;
        s = 3;
        BT = zeros(s+1,s+1);
        C = [1/2-sqrt(15)/10; 1/2; 1/2+sqrt(15)/10; 0];
        B = [0, 5/18, 4/9, 5/18];
        A = [5/36, 2/9-sqrt(15)/15, 5/36-sqrt(15)/30;
            5/36+sqrt(15)/24, 2/9, 5/36-sqrt(15)/24;
            5/36+sqrt(15)/30, 2/9+sqrt(15)/15, 5/36];
        BT(:,1) = C;
        BT(end,:) = B;
        BT(1:s,2:end) = A;
        order = 6;
        
    case 15
        solver = @odesImplicitRK;
        s = 2;
        BT = zeros(s+1,s+1);
        C = [1/2-sqrt(3)/6; 1/2+sqrt(3)/6; 0];
        B = [0, 1/2, 1/2];
        A = [1/4, 1/4-sqrt(3)/6;
            1/4+sqrt(3)/6, 1/4];
        BT(:,1) = C;
        BT(end,:) = B;
        BT(1:s,2:end) = A;
        order = 4;
        
    case 16
        solver = @odesImplicitRK;
        s = 5;
        BT = zeros(s+1,s+1);
        C = [0; (7-sqrt(21))/14; 1/2; (7+sqrt(21))/14; 1; 0];
        B = [0, 1/20, 49/180, 16/45, 49/180, 1/20];
        A = [0 0 0 0 0;
            1/14 1/9 (13-3*sqrt(21))/63 (14-3*sqrt(21))/126 0;
            1/32 (91+21*sqrt(21))/576 11/72 (91-21*sqrt(21))/576 0;
            1/14 (14+3*sqrt(21))/126 (13+3*sqrt(21))/63 1/9 0;
            0 7/18 2/9 7/18 0];
        BT(:,1) = C;
        BT(end,:) = B;
        BT(1:s,2:end) = A;
        order = 8;
        
    otherwise
        error('Invalid Solver Choice!');
        
end

%% Define the problem to be solved.

ODE = @(t,y) ex2(t,y);
yINI = [-1 2]'; % initial values
y_exact = @(x) exp(-x./2).*(-cos(2*x)+0.75*sin(2*x));
xMin = 0; xMax = 2;
x_gold = linspace(xMin,xMax,1000);

%% Convergence Analysis

N = 3:12;
t = length(N);
err_array = nan(1,t);
h_array = nan(1,t);
for s = 1:t
   Nx = 2^(N(s));
   x = linspace(xMin,xMax,Nx);
   h = x(2)-x(1);
   h_array(s) = h;
   y_ns = solver(ODE,xMin,xMax,Nx,yINI,BT);
   err_array(s) = norm(y_exact(x)-y_ns(1,:),2)./...
       norm(y_exact(x),2);
end

err_gage = h_array.^order;

figure(1)
loglog(h_array,err_array,'-b',...
    h_array,err_gage,'--r','linewidth',2);
title_str = sprintf('Convergence Test, Order = %d',order);
title(title_str,...
    'fontsize',14,'fontweight','bold');
xlabel('h','fontsize',12,'fontweight','bold');
ylabel('Relative Error','fontsize',12,'fontweight','bold');
grid on

%% Local Functions
function y = odesExplicitRK(ODE,a,b,N,yINI,BT)
% function y = odeExplicitRK(ODE,a,b,h,yINI,BT)
% y = solution (vector)
% ODE = function handle for y'
% a,b = begining and end of the interval for solution
% N = number of steps between a and b
% yINI = initial value for the solution
% BT = Butcher Tableau

% get Butcher Tableau Parameters
s = length(BT)-1;
c = BT(1:s,1);
B = BT(s+1,2:end);
A = BT(1:s,2:end);
stages = s;

x = linspace(a,b,N);
sys_size = length(yINI);
y = nan(sys_size,N);
y(:,1) = yINI;
h = x(2)-x(1);
for t = 1:(N-1)
    Xi = nan(sys_size,stages);
    
    for s = 1:stages
       Xi(:,s) = y(:,t);
       for i = 1:(s-1)
          Xi(:,s) = Xi(:,s) + h*A(s,i)*ODE(x(t)+c(i)*h,Xi(:,i)); 
       end
    end
    
    y(:,t+1) = y(:,t);
    for i = 1:stages
       y(:,t+1) = y(:,t+1) + h*B(i)*ODE(x(t)+c(i)*h,Xi(:,i)); 
    end
    
end
end

function y = odesImplicitRK(ODE,a,b,N,yINI,BT)

% get Butcher Tableau Parameters
s = length(BT) - 1;
c = BT(1:s,1);
b_t = BT(s+1,2:end);
A = BT(1:s,2:end);
[sys_size,~] = size(yINI);
x = linspace(a,b,N);
y = zeros(sys_size,N);
h = x(2) - x(1);
y(:,1) = yINI;

% SecantRootSys arguments
imax = 10000;
Err = 1e-16;
Ka = ones(sys_size,s)*.1; %<-- maybe zero is a bad choice...
Kb = ones(sys_size,s).*ODE(x(1),y(:,1));
for t = 1:(N-1)
    
   y(:,t+1) = y(:,t);
   
   % need to solve for the Ks
   F = @(iv,k) IRK_Ksol(ODE,iv,y(:,t),h,k,c,A);   
   K = SecantRootSys(F,x(t),Ka,Kb,imax,Err);
      
   for i = 1:s
       y(:,t+1) = y(:,t+1) + h*b_t(i)*ODE(x(t)+c(s)*h,K(:,i));
   end
   
   % update these for the secant solver for next time step
   Ka = Kb;
   Kb = K;
       
end
end

function R = IRK_Ksol(F,t,y,h,K_g,c,A)
% function K = IRK_Ksol(t,y,Kg,c,A) a generic function used to solve for K
% values of an IRK scheme
% Inputs
% F - function handle - of a function with two arguments (t,y) where t is
% the independent variable and y is a vector of dependent variables for the
% system
% t - scalar - current value of independent variable
% y - vector - current value of the dependent variables
% h - scalar - step size for independent variable
% K_g - matrix - guessed value of K - one column per k, one row per
% dependent variable
% c - vector - nodes for IRK scheme
% A - s x s matrix of coefficients for IRK scheme
%
% Output
% R  - residual matrix.  One column per k; one row per dependent variable

[s,~] = size(A); % number of stages for the IRK scheme
d = length(y); % number of dependent variables

R = nan(d,s);

for stage = 1:s
    R(:,stage) = y - K_g(:,stage);
    for i = 1:s
        R(:,stage) = R(:,stage) + h*(A(stage,i)*F(t+c(i)*h,K_g(:,i)));
    end    
end

% if K_g has all of the right K's for each stage and dependent variable,
% then R should be all zeros.
end

function Xs = SecantRootSys(Fun,t,Xa,Xb,imax,Err)
% function Xs = SecantRootSys(Fun,Xa,imax,Err) solves a system of
% non-linear equations using the Secant Method.
% 
% Inputs
% Fun - input function (must take 2 arguments: F(t,Xa) where t is a scalar
% representing the current value of the independent variable; and Xa is a vector containing the initial value for all dependent
% variables.
%
% 5 - scalar - initial value of independent variable
% Xa - vector - inital values
% Xb - vector - a second set of values
% imax - maximum number of iterations
% Err - error tolerance (relative 1-norm)

%Xb = Xa + Fun(t,Xa);

for i = 1:imax
   FunXb = Fun(t,Xb);
   % if Fun(t,x) is singular at Xb (i.e. Xb = 0), then I should
   % perturb Xb towards Xa
   if (isinf(FunXb))
      Xb = Xb - 0.1*(Xb-Xa);
      FunXb = Fun(t,Xb);
   end
   
   denom = Fun(t,Xa)-FunXb;
   
   % need to address problem of denom == 0
   % any zero elements of denom; change to 1
   denom(denom==0) = 1;
   
   Xi = Xb - FunXb.*(Xa - Xb)./(denom);
   
   % prevent further problems with zeros
   Xi(Xi==0) = rand;%<-- okay.  Really getting desparate now.
   
   
   if(norm(Xi - Xb,inf)/norm(Xb,inf)) < Err
       Xs = Xi;
       break;
   end
   
   Xa = Xb;
   Xb = Xi;
end

if i == imax
    error('Error!  Solution was not obtained at t=%g within %i iterations.',t,imax);
end
end

function [x,y] = odesEULER(ODEs,a,b,N,yINI,~)
[sys_size,~] = size(yINI);
h = (b - a)/N; % interval size
x = linspace(a,b,N+1);
y = zeros(sys_size,N+1);

y(:,1) = yINI;
for i = 1:N
    y(:,i+1) = y(:,i) + ODEs(x(i),y(:,i))*h;
end
end

function dw = ex2(~,w) % generally expect 2 arguments for solvers (IV first)
dw = nan(2,1);
dw(1) = w(2);
dw(2) = -0.25*(4*w(2) + 17*w(1));
end
