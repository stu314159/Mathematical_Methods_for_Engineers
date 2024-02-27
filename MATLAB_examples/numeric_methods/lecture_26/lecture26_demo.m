%% Lecture 26 - MATLAB Demo
% a comprehensive demonstration of RK methods
clear
clc
close 'all'

%% Pick an ODE
ODE_choice = 2;


switch ODE_choice
    % A collection of 2nd order IVPs.
    case 1
        ODE = @ex10_7;
        a = 0; b = 3;
        yINI = [3; 0.2];
                
    case 2 % nice oscillatory behavior
        ODE = @ex10_8;
        a = 0; b = 18;
        yINI = [pi/2; 0];
        
    case 3
        ODE = @Problem10p30;
        a = 0; b = 0.1;
        yINI = [0;0];
        
    case 4
        ODE = @Problem10p6;
        a = 0; b = 1.2;
        yINI = [1; 1];
        
    case 5
        ODE = @stiffODE;
        a = 0;
        b = 0.5;
        yINI = [1;2];
        
end

%% Set Discretization
N = 500;
x = linspace(a,b,N);


%% Pick a Solver
Solver_choice = 9;

% Explicit Methods
% 1 = 1st order, explicit Euler
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
% 10 = 5th order, 3-stage implicit RK "Radeau IIA" 
% 11 = 4th order Lobatto IIIC (try for stiff problems) 
% 12 = 1 stage Backward Euler 
% 13 = 2nd-order, implicit midpoint 
% 14 = Gauss Legendre method of order 6 (looks good)
% 15 = 4th order IRK from Butcher, Eq 237c (looks good)

switch Solver_choice
    
    case 1
        solver = @odesEULER;
        BT = [];
        
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
        
    case 7
        solver = @odesExplicitRK;
        s = 6;
        BT = zeros(s+1,s+1);
        C = [0; 1/4; 3/8; 12/13; 1; 1/2; 0];
        B = [0 16/135 0 6656/12825 28561/56430 -9/50 2/55];
        %B = [0 25/216 0 1408/2565 2197/4104 -1/5 0];
        A = [0 0 0 0 0 0;
            1/4 0 0 0 0 0;
            3/32 9/32 0 0 0 0;
            1932/2197 -7200/2197 7296/2197 0 0 0;
            439/216 -8 3680/513 -845/4104 0 0;
            -8/27 2 -3544/2565 1859/4104 -11/40 0;];
        BT(:,1) = C;
        BT(end,:) = B;
        BT(1:s,2:end) = A;
        
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
        
    case 10
        % 4-stage implicit RK "Radeau IIA" 
        solver = @odesImplicitRK;
        s = 3;
        BT = zeros(s+1,s+1);
        C = [2/5 - sqrt(6)/10; 2/5+sqrt(6)/10; 1; 0];
        B = [0, 4/9 - sqrt(6)/36, 4/9 + sqrt(6)/36, 1/9];
        A = [11/45 - 7*sqrt(6)/360, 37/225 - 169*sqrt(6)/1800, -2/225+sqrt(6)/75;
            37/225+169*sqrt(6)/1800, 11/45+7*sqrt(6)/360, -2/225-sqrt(6)/75;
            4/9 - sqrt(6)/36, 4/9+sqrt(6)/36, 1/9];
        BT(:,1) = C;
        BT(end,:) = B;
        BT(1:s,2:end) = A;
        
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


%% Compute Numeric Solution

y = solver(ODE,a,b,N,yINI,BT);

%% plot a solution
figure(1)
plot(x,y(1,:),'-r',...
    x,y(2,:),'-g',...
    'linewidth',3);
%plot(x,y,'-r','linewidth',3);
grid on;
title('Numeric Solution','fontsize',16,'fontweight','bold');
xlabel('X','fontsize',14,'fontweight','bold');
ylabel('Y','fontsize',14,'fontweight','bold');
set(gca,'fontsize',12,'fontweight','bold');

%% Compare with MATLAB Built-in Method

M_RT = 1e-13;  M_AT = 1e-11;
M = ode; M.ODEFcn = ODE; M.InitialValue = yINI;
M.InitialTime = a;
M.RelativeTolerance = M_RT;
M.AbsoluteTolerance = M_AT;

refSol = solve(M,x);

figure(2)
plot(refSol.Time,refSol.Solution(1,:),...
    '-r',...
    refSol.Time,refSol.Solution(2,:),'-g',...
    'linewidth',3,'markersize',4);
title('Solution with ODE','fontsize',16,'fontweight','bold');
grid on
xlabel('T','fontsize',14,'fontweight','bold');
ylabel('Y(T)','fontsize',14,'fontweight','bold');
set(gca,'fontsize',12,'fontweight','bold');

%% Compare Selected Method with MATLAB Method
M1_interp = interp1(refSol.Time,refSol.Solution(1,:),x);
M2_interp = interp1(refSol.Time,refSol.Solution(2,:),x);

rel_err1 = norm(M1_interp - y(1,:),2)./norm(M1_interp,2);
rel_err2 = norm(M2_interp - y(2,:),2)./norm(M2_interp,2);

fprintf('Relative error for y1: %g \n',rel_err1);
fprintf('Relative error for y2: %g \n',rel_err2);

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

function dF = Problem10p6(t,x)
% x is a vector; x(1) corresponds to "x", x(2) corresponds to "y"
dF = nan(2,1);
dF(1) = x(1) - x(2).*t;
dF(2) = t + x(2);

end

function df = ex10_7(x,y)
df = zeros(size(y));
df(1) = (-y(1) + y(2))*exp(1-x) + 0.5*y(1);
df(2) = y(1) - y(2).^2;


end

function df = ex10_8(~,y)
c = 0.16; % N-s/m, damping coefficient
L = 1.2; % m, length
m = 0.5; % kg, mass
g = 9.81; % m/s^2, acceleration due to graviation

df = zeros(size(y));
df(1) = y(2);
df(2) = -(c/m)*y(2) - (g/L)*sin(y(1));


end

function dF = Problem10p30(x,y)
% problem-specific constants
V = 165;
N = 600;
omega = 120*pi;
a = 0.14;
omega_o = 83;

dF = nan(2,1);
dF(1) = y(2);
dF(2) = -(omega*V/N)*cos(omega*x) - (omega_o^2)*y(1) - a*(y(1).^3);

end

function dF = stiffODE(~,x)
dF = nan(2,1);
dF(1) = 998*x(1) - 1998*x(2);
dF(2) = 1000*x(1) - 2000*x(2);
end

function y = odesImplicitRK(ODE,a,b,N,yINI,BT)

% get Butcher Tableau Parameters
s = length(BT) - 1;
c = BT(1:s,1);
b_t = BT(s+1,2:end);
A = BT(1:s,2:end);
[sys_size,~] = size(yINI);
%h = (b-a)/N;
x = linspace(a,b,N);
h = x(2) - x(1);
y = zeros(sys_size,N);
y(:,1) = yINI;

% SecantRootSys arguments
imax = 10000;
Err = 1e-14;
Ka = ones(sys_size,s)*.01; %<-- maybe zero is a bad choice...
Kb = ones(sys_size,s).*ODE(x(1),y(:,1));
for t = 1:(N-1)
    
    y(:,t+1) = y(:,t);

    % need to solve for the Ks
    F = @(iv,k) IRK_Ksol(ODE,iv,y(:,t),h,k,c,A);
    K = SecantRootSys(F,x(t),Ka,Kb,imax,Err);
    % options = optimoptions('fsolve',...
    %     'Display','None',...
    %     'FunctionTolerance',1e-10);
    % K = fsolve(F,Kb,options);

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
    R(:,stage) = y; 
    for i = 1:s
        R(:,stage) = R(:,stage) + ...
            h*(A(stage,i)*F(t+c(i)*h,K_g(:,i)));
    end
    R(:,stage) = R(:,stage) - K_g(:,stage);
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
% t - scalar - initial value of independent variable
% Xa - vector - inital values
% Xb - vector - a second set of values
% imax - maximum number of iterations
% Err - error tolerance (relative 1-norm)

%Xb = Xa + Fun(t,Xa);

for i = 1:imax
   FunXb = Fun(t,Xb);
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

function y = odesEULER(ODEs,a,b,N,yINI,~)
[sys_size,~] = size(yINI);
h = (b - a)/N; % interval size
x = linspace(a,b,N);
y = zeros(sys_size,N);

y(:,1) = yINI;
for i = 1:(N-1)
    y(:,i+1) = y(:,i) + ODEs(x(i),y(:,i))*h;
end
end