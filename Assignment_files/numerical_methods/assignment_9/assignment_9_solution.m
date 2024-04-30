%% Assignment #9 Solution
clear
clc
close 'all'

%% Problem #1
fprintf('\n\n Problem #1 \n\n');
a = 0; b = pi;
Fx = @(x) 2*x;
Gx = @(x) 5;
Hx = @(x) cos(3*x);
Ya = 1.5; Yb = 0;
n = 100;
Wa = -5; Wb = -1.5;
[x,y] = BVPShootSecant(Fx,Gx,Hx,a,b,n,Ya,Yb,Wa,Wb);

%% Plot the solution
figure(1)
plot(x,y(1,:),'-b','linewidth',3);
title('Problem #1 Solution','fontsize',16,'fontweight','bold')
xlabel('X','fontsize',14,'fontweight','bold');
ylabel('Y(X)','fontsize',14,'fontweight','bold');
grid on
set(gca,'fontsize',12,'fontweight','bold');


%% Problem #2
fprintf('\n\n Problem #2 \n\n');
fprintf('Part a) \n');
a = 1; b = 3.5; n = 100;
Fx = @(x) 1./x;
Gx = @(x) 0;
Hx = @(x) -500./x;
Ya = 600; Yb = 25;
Wa = -1; Wb = 0;

[x2,y2] = BVPShootSecant(Fx,Gx,Hx,a,b,n,Ya,Yb,Wa,Wb);

%% Plot the solution
figure(2)
plot(x2,y2(1,:),'-b','linewidth',3);
title('Problem #2.a BVPShootSecant','fontsize',16,'fontweight','bold')
xlabel('X','fontsize',14,'fontweight','bold');
ylabel('Y(X)','fontsize',14,'fontweight','bold');
grid on
set(gca,'fontsize',12,'fontweight','bold');

fprintf('Slope at r = %g using BVPShootSecant: %g \n',...
    x2(end),y2(2,end));

%% Use BVP5C
ode2 = @(r,T) [T(2); (-1./r).*T(2) - 500./r];
bcfun = @(ya,yb) [ya(1) - Ya; yb(1) - Yb];
Yguess = [Ya 0];
solinit = bvpinit(linspace(a,b,n),Yguess);
sol2a = bvp5c(ode2,bcfun,solinit);

%% Plot the solution
figure(3)
plot(sol2a.x,sol2a.y(1,:),'-b','linewidth',3);
title('Problem #2.a BVP4C','fontsize',16,'fontweight','bold')
xlabel('X','fontsize',14,'fontweight','bold');
ylabel('Y(X)','fontsize',14,'fontweight','bold');
grid on
set(gca,'fontsize',12,'fontweight','bold');

fprintf('Slope at r = %g using BVP4C: %g \n',...
    sol2a.x(end),sol2a.y(2,end));

%% Problem #2, part b - experiment with RHS
fprintf('Part b) \n');
RHS = -1500;
Hx = @(x) RHS./x;
[x2b,y2b] = BVPShootSecant(Fx,Gx,Hx,a,b,n,Ya,Yb,Wa,Wb);

figure(4)
plot(x2b,y2b(1,:),'-b','linewidth',3);
title('Problem #2.b BVPShootSecant','fontsize',16,'fontweight','bold')
xlabel('X','fontsize',14,'fontweight','bold');
ylabel('Y(X)','fontsize',14,'fontweight','bold');
grid on
set(gca,'fontsize',12,'fontweight','bold');


%% Problem #2.c change boundary conditions
fprintf('Part c) \n');
T_inf = 25; % C, ambient temp
h = 1.5; % W/cm^2-K

bcfun = @(ya,yb) [ya(1) - Ya; yb(2) + h.*(yb(1)-T_inf)];
Yguess = [Ya 0];
solinit = bvpinit(linspace(a,b,n),Yguess);
sol2c = bvp5c(ode2,bcfun,solinit);

figure(5)
plot(sol2c.x,sol2c.y(1,:),'-b','linewidth',3);
title('Problem #2.c BVP4C - Mixed BC','fontsize',16,'fontweight','bold')
xlabel('X','fontsize',14,'fontweight','bold');
ylabel('Y(X)','fontsize',14,'fontweight','bold');
grid on
set(gca,'fontsize',12,'fontweight','bold');

fprintf('T at r = %g: %g C \n',sol2c.x(end),sol2c.y(1,end));

%% Problem #3
fprintf('\nProblem #3 \n\n');

rho_e = 32e-8; % Ohm-meter, electrical resistivity
I = 0.5; % A, current
D = 1e-4; % m, 
k = 72; % W/m/K, thermal conductivity
C = -I^2*rho_e/(k*(0.25*pi*D^2)^2);
rMin = 1e-6;
rMax = D/2;
n = 50;

ode3 = @(r,T) [T(2); (-1./r).*T(2) + C];
bcfun_a = @(Ta,Tb) [Ta(2)-0; Tb(1)-300];
Tguess = [500 0];
solinit = bvpinit(linspace(rMin,rMax,n),Tguess);

sol3a = bvp5c(ode3,bcfun_a,solinit);

%% Plot the solution for Part A

figure(6)
plot(sol3a.x,sol3a.y(1,:),'-b','linewidth',3)
title('Temperature Profile Prob 3.a',...
    'fontsize',14,'fontweight','bold');
xlabel('R [m]','fontsize',12,'fontweight','bold');
ylabel('T(R) [K]','fontsize',12,'fontweight','bold');
grid on
set(gca,'fontsize',10,'fontweight','bold');

fprintf('Temperature difference 3.a = %g K \n',...
    sol3a.y(1,1) - sol3a.y(1,end));
%%
h = 100; % W/m^2-K, convective heat transfer coefficient
T_inf = 300; % K, ambient temperature
bcfun_b = @(Ta,Tb) [Ta(2)-0; Tb(2) + (h/k)*(Tb(1) - T_inf)]; 

sol3b = bvp5c(ode3,bcfun_b,solinit);

figure(7)
plot(sol3b.x,sol3b.y(1,:),'-b','linewidth',3)
title('Temperature Profile Prob 3.b',...
    'fontsize',14,'fontweight','bold');
xlabel('R [m]','fontsize',12,'fontweight','bold');
ylabel('T(R) [K]','fontsize',12,'fontweight','bold');
grid on
set(gca,'fontsize',10,'fontweight','bold');

fprintf('Temperature difference 3.b = %g K \n',...
    sol3b.y(1,1) - sol3b.y(1,end));
%% Local Functions
function [x,y] = BVPShootSecant(Fx,Gx,Hx,a,b,n,Ya,Yb,Wa,Wb)
% function [x,y] = BVPShootSecant(Fx,Gx,Hx,a,b,n,Ya,Yb,Wa,Wb)
% solves a second-order BVP of form:
% y'' + f(x)y' + g(x)y = h(x) on the domain a < x < b with 
% y(a) = Ya, and y(b) = Yb.
%
% inputs
% Fx - function handle
% Gx - function handle
% Hx - function handle
% a - scalar domain left boundary
% b - scalar domain right boundary
% n - scalar, number of subdivisions
% Ya - scalar, left BC
% Yb - scalar, right BC
% Wa - scalar slope low guess
% Wb - scaalr slope high guess for bisection

% define ODE: x - independent variable, y - vector of 2 dependent
% variables.
ODE = @(x,y) [y(2); -Fx(x).*y(2) - Gx(x).*y(1) + Hx(x)];

[~,yL] = odesCRK4(ODE,a,b,n,[Ya; Wa]);
[~,yH] = odesCRK4(ODE,a,b,n,[Ya; Wb]);


Ei = yL(1,end) - Yb;
Eim1 = yH(1,end) - Yb;
Wim1 = Wb;
Wi = Wa;



%% set iteration parameters
tol = 1e-3;
imax = 100;

for i = 1:imax
   Wip1 = Wi - (Wim1 - Wi)*Ei/(Eim1 - Ei);
   [x,y] = odesCRK4(ODE,a,b,n,[Ya; Wip1]);
   Eim1 = Ei;
   Ei = y(1,end) - Yb;
   
   if abs(Ei) < tol
      break; %<-- success! 
   end
   
   Wim1 = Wi;
   Wi = Wip1;
end

if i == imax
    fprintf('Failed to find a solution within error tolerance in %d iterations\n',...
        imax);
end
end

function [t,y] = odesCRK4(F,a,b,N,yINI)
assert(min(size(yINI))==1);

% from rk demo script
stages = 4;
c = [0; 1/2; 1/2; 1];
B = [1/6 2/6 2/6 1/6];
A = [0 0 0 0;
    1/2 0 0 0;
    0 1/2 0 0;
    0 0 1 0;];


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
          Xi(:,s) = Xi(:,s) + h*A(s,i)*F(x(t)+c(i)*h,Xi(:,i)); 
       end
    end
    
    y(:,t+1) = y(:,t);
    for i = 1:stages
       y(:,t+1) = y(:,t+1) + h*B(i)*F(x(t)+c(i)*h,Xi(:,i)); 
    end
    
end

t = x;
end