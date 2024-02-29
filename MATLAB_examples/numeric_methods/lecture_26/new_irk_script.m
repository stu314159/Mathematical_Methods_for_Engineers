%% New RK
clear
clc
close 'all'

%% Select Solver
Solver_choice = 10;

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

%% Define the problem to be solved
% start with an autonomous problem.

F = @(t,y) y; % dy/dt = y
yINI = 1; Nt = 5;
y_exact = @(t) exp(t);
tMin = 0; tMax = 5;
tGold = linspace(tMin,tMax,1000);
t = linspace(tMin,tMax,Nt);

%% Solve
y = odesImplicitRK(F,tMin,tMax,Nt,yINI,BT);

%% Plot and Compare
figure(1)
plot(tGold,y_exact(tGold),'-r',...
    t,y,'sb',...
    'markersize',10,'linewidth',3);
grid on
title('Sovler Test')


%% Local function
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

for ts = 1:(N-1)
    
   y(:,ts+1) = y(:,ts);
   
   % need to solve for the Xis
   Xi = find_xi(ODE,x(ts),y(:,ts+1),h,A,c);
      
   for i = 1:s
       y(:,ts+1) = y(:,ts+1) + h*b_t(i)*ODE(x(ts)+c(s)*h,Xi(:,i));
   end
   
end

end

function xi = find_xi(ODE,t,y,h,A,c)

[s,~] = size(A); % get number of stages
[n] = length(y); % get number of dependent variables
options = optimoptions('fsolve','Display','none');
xi = nan(n,s); % pre-allocate the array for xi
for v = 1:n % iterate over the dependent variables
    xi_end_est = y+ODE(t,y)*h; % if slope constant over h, what is last xi
    xi_0 = linspace(y,xi_end_est,s); % initial estimate of sample points
    R = @(k) y+h*(A*ODE(t+c*h,k')) - k';
    xi_n = fsolve(R,xi_0,options);
    xi(v,:) = xi_n';
end

end

function dw = ex2(~,w) % generally expect 2 arguments for solvers (IV first)
dw = nan(2,1);
dw(1) = w(2);
dw(2) = -0.25*(4*w(2) + 17*w(1));
end