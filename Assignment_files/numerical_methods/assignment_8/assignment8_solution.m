%% Assignment #8 Solution
clear
clc
close 'all'

%% Problem #1
fprintf('\n Problem #1 \n\n');
f = @(x,y,z) x.^2 - y.^2 - z.^2;

int3_sol = integral3(f,-1,1,-1,1,-1,1);
fprintf('Solution with integral3: %16.15f \n',int3_sol);

% trial area
xMin = -1; xMax = 1;
yMin = -1; yMax = 1;
zMin = -1; zMax = 1;
Vol = (xMax-xMin)*(yMax-yMin)*(zMax-zMin);

s = 10:26; % no need to go overboard on the numbers
num_steps = length(s);
mc_est = nan(1,num_steps);
N = 2.^s;

for i = 1:num_steps
    fprintf('Computing with %g samples.\n',N(i));
    x = xMin + (xMax-xMin)*rand(N(i),1);
    y = yMin + (yMax-yMin)*rand(N(i),1);
    z = zMin + (zMax-zMin)*rand(N(i),1);
    mc_est(i) = (Vol/N(i))*sum(f(x,y,z)); 
end

rel_err = abs(int3_sol - mc_est)./abs(int3_sol);
err_gage = sqrt(1./N);
%% plot the results
figure(1)
loglog(N,rel_err,'-b',...,
    N,err_gage,'--r','linewidth',3);
grid on
title('Monte Carlo Convergence',...
    'fontsize',14,'fontweight','bold');
xlabel('Number of Samples',...
    'fontsize',12,'fontweight','bold');
ylabel('Relative Error',...
    'fontsize',12,'fontweight','bold');
set(gca,'fontsize',10,'fontweight','bold');
legend('Monte Carlo Error','N^{-1/2} Convergence Gage');

%% Problem #2
fprintf('\n Problem #2 \n\n');

% Major axis dimensions
a = 9.5/2; % cm
b = 8/2; % cm
c = 4.2/2; % cm

delta = 1 - (c/a)^2;
epsilon = 1 - (c/b)^2;
p = @(phi) delta*sin(phi).^2 + ...
    epsilon*cos(phi).^2;
F = @(theta,phi) sin(theta).*...
    sqrt(1-p(phi).*sin(theta).^2);

thetaMin = 0; thetaMax = pi/2;
phiMin=0; phiMax = pi/2;
p2_int = 8*a*b*integral2(F,...
    thetaMin,thetaMax,phiMin,phiMax);

fprintf('Tumor surface area = %g cm^3. \n',...
    p2_int);

%% Problem #3
fprintf('\n Problem #3 \n\n');

F = @(x,y) prob3(x,y);
y_exact = @(x) exp(-x).*(-cos(x)-4*sin(x)./5);
a = 0; b = 1.5;
yINI=[-1; 0.2];
N1 = 100;
[x1,y1] = odesEULER(F,a,b,N1,yINI);
rel_err1 = norm(y_exact(x1)-y1(1,:),2)/norm(y_exact(x1),2);
N2 = 1000;
[x2,y2] = odesEULER(F,a,b,N2,yINI);
rel_err2 = norm(y_exact(x2)-y2(1,:),2)/norm(y_exact(x2),2);

fprintf('N = %d, relative error: %g \n',N1,rel_err1);
fprintf('N = %d, relative error: %g \n',N2,rel_err2);
fprintf('N2/N1 = %g \n',N2/N1);
fprintf('Error 1 / Error 2 = %g \n',rel_err1/rel_err2);

%% Problem #4 (MATLAB check)
fprintf('\n Problem #4 \n\n');
F = @(t,y) y + t.^3;
yINI = 1;
tMin = 0;
tMax = 1.5;
N = 4;
[t4,y4] = odesCRK4(F,tMin,tMax,N,yINI);
fprintf('y(%g) = %g \n',t4(2),y4(2));
fprintf('y(%g) = %g \n',t4(3),y4(3));
fprintf('y(%g) = %g \n',t4(4),y4(4));

y4_ex = @(t) 7*exp(t) - t.^3 - 3*t.^2 - 6*t - 6;

err4 = y4_ex(t4) - y4;
fprintf('err(%g) = %g \n',t4(2),err4(2));
fprintf('err(%g) = %g \n',t4(3),err4(3));
fprintf('err(%g) = %g \n',t4(4),err4(4));

    
%% Problem #5
fprintf('\n Problem #5 \n\n');
A_tank = 3.13; % m^2, cross sectional area of tank
%A_tank = 3.13; % adjust this to play
A_pipe = 0.06; % m^2, cross sectional area of drain pipe
%A_pipe = 0.1; % adjust this to play
C = pi/12;
K1 = 300; % kg/s
%K1 = 600; % adjust this to play
%K2 = 1000; % kg/s
K2 = 2000; % adjust this to play
rho = 1000; % kg/m^3
g = 9.81; % m/s^2

F = @(t,h) (1./(rho*A_tank))*...
    (K1+K2*sin(5*C*t).*cos(C*t) - ...
    rho*A_pipe*sqrt(2*g*h));

tMin = 0; tMax = 150; % seconds
N = 1000;
hINI = 3; % m, initial height of water in the tank.

[t_Euler,y_Euler] = odesEULER(F,tMin,tMax,N,hINI);

% Plot the results
figure(3)
plot(t_Euler,y_Euler,'-b','linewidth',3);
title('Problem #5 - Euler Method',...
    'fontsize',14,'fontweight','bold');
xlabel('time [sec]','fontsize',12,'fontweight','bold');
ylabel('height [m]','fontsize',12,'fontweight','bold');
grid on
set(gca,'fontsize',10,'fontweight','bold');

% output requested value
fprintf('(Euler) height in tank at t=%g seconds is %g m.\n',...
    t_Euler(end),y_Euler(end));

%% Solve again, but using odesCRK4
[t_CRK4,y_CRK4] = odesCRK4(F,tMin,tMax,N,hINI);

% Plot the results
figure(4)
plot(t_CRK4,y_CRK4,'-b','linewidth',3);
title('Problem #5 - Classical RK4 Method',...
    'fontsize',14,'fontweight','bold');
xlabel('time [sec]','fontsize',12,'fontweight','bold');
ylabel('height [m]','fontsize',12,'fontweight','bold');
grid on
set(gca,'fontsize',10,'fontweight','bold');

% output requested value
fprintf('(CRK4) height in tank at t=%g seconds is %g m.\n',...
    t_CRK4(end),y_CRK4(end));

%% Problem #6
fprintf('\n Problem #6 \n\n');
F = @(t,y) prob6(t,y);
N = 1000;
yINI = [0; 0];
tMin = 0; tMax = 3;

[t_6,y_6] = odesCRK4(F,tMin,tMax,N,yINI);
dist = y_6(1,:);
vel = y_6(2,:);
acc = FirstDeriv(t_6,vel);

figure(5)

subplot(3,1,1)

plot(t_6,dist,'-c','linewidth',3);
title('Problem #6 Solution','fontsize',14,...
    'fontweight','bold');
ylabel('distance [ft]','fontsize',12,'fontweight','bold');
set(gca,'fontsize',10,'fontweight','bold');
grid on;

subplot(3,1,2)
plot(t_6,vel,'-c','linewidth',3);
ylabel('velocity [ft/s]','fontsize',12,'fontweight','bold');
set(gca,'fontsize',10,'fontweight','bold');
grid on

subplot(3,1,3)
plot(t_6,acc,'-c','linewidth',3);
ylabel('acceleration [ft/s^2]','fontsize',12,'fontweight','bold');
xlabel('time [s]','fontsize',12,'fontweight','bold');
set(gca,'fontsize',10,'fontweight','bold');
grid on

fprintf('At t = 3 seconds:\n');
fprintf('distance = %g ft \n',dist(end));
fprintf('velocity = %g ft/s \n',vel(end));
fprintf('acceleration = %g ft/s \n',acc(end));

snapnow


%% Local Functions
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

function [t,y] = odesEULER(F,a,b,N,yINI)

assert(min(size(yINI))==1);
ndof = length(yINI);
y = nan(ndof,N); % construct output array
t = linspace(a,b,N);
h = t(2)-t(1);
y(:,1) = yINI;

for i = 1:(N-1)
   y(:,i+1) = y(:,i)+F(t(i),y(:,i))*h; 
end

end

function dw = prob3(~,w)
dw = nan(2,1);
dw(1) = w(2);
dw(2) = -2*w(2) - 2*w(1);
end

function dz = prob6(t,z)
g = 32.2; % ft/s^2, gravitational acceleration
%D = @(t) 0.005*g*z(2).^2;% drag
D = @(t) 0.0025*g*z(2).^2; % play with this.
w = @(t) 3000 - 80*t;
%w = @(t) 3000; % play with this
T = 8000; % lb
%T = 20000; % play with this.

dz = nan(2,1);
dz(1) = z(2);
dz(2) = (g./w(t))*T - g - D(t)*g./w(t);

end

function yd = FirstDeriv(x,y)
assert(min(size(x))==1);% x and y should both be vectors
assert(min(size(y))==1);
assert(max(size(x))==max(size(y))); % x and y the same length
assert(max(size(y))>=3); % enough data for proposed method
[m,n] = size(y);
yd = nan(m,n); % construct output same shape as y.
h = x(2)-x(1); % it would be good to verify this
for s = 2:(length(h)-1)
    ht = x(s+1)-x(s);
    assert(abs(h - ht)<=eps);
end
% but might be going overboard.

% left end
yd(1) = (-3*y(1) + 4*y(2) - y(3))/(2*h);

% right end
yd(end) = (y(end-2) - 4*y(end-1) + 3*y(end))/(2*h);

% middle
i = 2:(length(x)-1);
yd(i) = (y(i+1) - y(i-1))/(2*h);

end

% function y = odesExplicitRK(ODE,a,b,N,yINI,BT)
% function y = odeExplicitRK(ODE,a,b,h,yINI,BT)
% y = solution (vector)
% ODE = function handle for y'
% a,b = begining and end of the interval for solution
% N = number of steps between a and b
% yINI = initial value for the solution
% BT = Butcher Tableau
% 
% get Butcher Tableau Parameters
% s = length(BT)-1;
% c = BT(1:s,1);
% B = BT(s+1,2:end);
% A = BT(1:s,2:end);
% stages = s;
% 
% x = linspace(a,b,N);
% sys_size = length(yINI);
% y = nan(sys_size,N);
% y(:,1) = yINI;
% h = x(2)-x(1);
% for t = 1:(N-1)
%     Xi = nan(sys_size,stages);
%     
%     for s = 1:stages
%        Xi(:,s) = y(:,t);
%        for i = 1:(s-1)
%           Xi(:,s) = Xi(:,s) + h*A(s,i)*ODE(x(t)+c(i)*h,Xi(:,i)); 
%        end
%     end
%     
%     y(:,t+1) = y(:,t);
%     for i = 1:stages
%        y(:,t+1) = y(:,t+1) + h*B(i)*ODE(x(t)+c(i)*h,Xi(:,i)); 
%     end
%     
% end
% end
% 
% 
