%% Lecture 23 Demo
clear
clc
close 'all'

%% Define the problem to be solved.
%%
% 
% $$y^{\prime} = x^2/y, \ \ y(0)=2$$
% 
f = @(x,y) (x.^2)./y;
y_exact = @(x) sqrt((2/3)*x.^3 + 4);
xMin = 0; xMax = 2.0;


%% Euler Explicit Demonstration

N = 30;
x = linspace(xMin,xMax,N+1);
y_ns = nan(1,N);
y_ns(1) = 2;
h = x(2)-x(1);
for t = 1:N
   y_ns(t+1) = y_ns(t)+f(x(t),y_ns(t))*h; 
end

%% Plot the Solution
x_gold = linspace(xMin,xMax,1000);

figure(1)
plot(x,y_ns,'-.b',...
    x_gold,y_exact(x_gold),'-r',...
    'linewidth',3);
title("Solution with Euler's Method",'fontsize',14,...
    'fontweight','bold');
xlabel('X','fontsize',12,'fontweight','bold');
ylabel('Y(X)','FontSize',12,'FontWeight','bold');
grid on;
set(gca,'fontsize',10,'fontweight','bold');
legend('Euler Explicit','Exact Solution');

%% Euler Explicit Convergence Analysis

N = 3:18;
t = length(N);
err_array = nan(1,t);
h_array = nan(1,t);
for s = 1:t
   Nx = 2^(N(s));
   x = linspace(xMin, xMax, Nx);
   h = x(2)-x(1);
   h_array(s) = h;
   y_ns = nan(1,Nx); y_ns(1) = 2;
   for i = 1:(Nx-1)
       y_ns(i+1)=y_ns(i)+f(x(i),y_ns(i))*h;
   end
   err_array(s) = norm(y_exact(x)-y_ns,2)./...
       norm(y_exact(x),2);
end

err_gage = h_array;

figure(2)
loglog(h_array,err_array,'-b',...
    h_array,err_gage,'--r','linewidth',2);
title("Convergence Euler's Explicit Method",...
    'fontsize',14,'fontweight','bold');
xlabel('h','fontsize',12,'fontweight','bold');
ylabel('Relative Error','fontsize',12,'fontweight','bold');
grid on
legend('Convergence Rate','h');
%% Forward Euler's Method

N = 30;
x = linspace(xMin,xMax,N);
y_ns = nan(1,N);
y_ns(1) = 2;
h = x(2)-x(1);
options = optimoptions('fsolve','Display','none');

for t = 1:N-1
   fe_fun = @(y) y - y_ns(t) - f(x(t+1),y)*h;   
   y_ns(t+1) = fsolve(fe_fun,y_ns(t),options);
end

figure(3)
plot(x,y_ns,'-.b',...
    x_gold,y_exact(x_gold),'-r',...
    'linewidth',3);
title("Solution with Forward Euler's Method",'fontsize',14,...
    'fontweight','bold');
xlabel('X','fontsize',12,'fontweight','bold');
grid on;
set(gca,'fontsize',10,'fontweight','bold');

%% A Stiff Example with Euler's (explicit) Method
% an example Stiff equation from Trefethen
% Finite Difference and Spectral Methods section 1.8
f_stiff = @(x,y) -100*(y-cos(x))-sin(x);
f_stiff_exact = @(x,y) cos(x);

N = 95; % see what happens if N is made smaller.
x = linspace(xMin,xMax,N+1);
y_ns = nan(1,N+1);
y_ns(1) = 1;
h = x(2)-x(1);
for t = 1:N
   y_ns(t+1) = y_ns(t)+f_stiff(x(t),y_ns(t))*h; 
end

figure(4)
plot(x,y_ns,'-.b',...
    x_gold,f_stiff_exact(x_gold),'--r',...
    'linewidth',3);
title("Stiff Equation with Euler's Method",'fontsize',14,...
    'fontweight','bold');
xlabel('X','fontsize',12,'fontweight','bold');
ylabel('Y(X)','FontSize',12,'FontWeight','bold');
grid on;
set(gca,'fontsize',10,'fontweight','bold');
legend('Euler Explicit','Exact Solution');

%% A Stiff Example with Euler's (Implicit) Method
N = 10;
x = linspace(xMin,xMax,N+1);
y_ns = nan(1,N);
y_ns(1) = 1;
h = x(2)-x(1);
options = optimoptions('fsolve','Display','none');
for t = 1:N
   fe_fun = @(y) y - y_ns(t) - f_stiff(x(t+1),y)*h;
   y_ns(t+1) = fsolve(fe_fun,y_ns(t),options);
end

figure(5)
plot(x,y_ns,'-.b',...
    x_gold,f_stiff_exact(x_gold),'--r',...
    'linewidth',3);
title("Stiff Equation with Euler's Method",'fontsize',14,...
    'fontweight','bold');
xlabel('X','fontsize',12,'fontweight','bold');
grid on;
set(gca,'fontsize',10,'fontweight','bold');
ylabel('Y(X)','FontSize',12,'FontWeight','bold');
grid on;
set(gca,'fontsize',10,'fontweight','bold');
legend('Euler Implicit','Exact Solution');


%% Modified Euler's Method

N = 30;
x = linspace(xMin,xMax,N+1);
y_ns = nan(1,N+1);
y_ns(1) = 2;
h = x(2)-x(1);

for t = 1:N
   f_xy = f(x(t),y_ns(t));
   % initial Euler step
   y_ns_EU = y_ns(t) + f_xy*h;
   
   % Apply Mod Euler
   y_ns(t+1) = y_ns(t) + ...
       (f_xy + f(x(t+1),y_ns_EU))*h/2;   
end

figure(6)
plot(x,y_ns,'-.b',...
    x_gold,y_exact(x_gold),'--r',...
    'linewidth',3);
title("Solution with Modified Euler's Method",'fontsize',14,...
    'fontweight','bold');
xlabel('X','fontsize',12,'fontweight','bold');
grid on;
set(gca,'fontsize',10,'fontweight','bold');

%% Modified Euler Convergence Analysis

N = 3:18;
t = length(N);
err_array = nan(1,t);
h_array = nan(1,t);
for s = 1:t
   Nx = 2^(N(s));
   x = linspace(xMin, xMax, Nx);
   h = x(2)-x(1);
   h_array(s) = h;
   y_ns = nan(1,Nx); y_ns(1) = 2;
   for i = 1:(Nx-1)
       f_xy = f(x(i),y_ns(i));
       % initial Euler step
       y_ns_EU = y_ns(i) + f_xy*h;
       
       % Apply Mod Euler
       y_ns(i+1) = y_ns(i) + ...
           (f_xy + f(x(i+1),y_ns_EU))*h/2;
   end
   err_array(s) = norm(y_exact(x)-y_ns,2)./...
       norm(y_exact(x),2);
end

err_gage = h_array.^2;

figure(7)
loglog(h_array,err_array,'-b',...
    h_array,err_gage,'--r','linewidth',2);
title("Convergence Modified Euler's Method",...
    'fontsize',14,'fontweight','bold');
xlabel('h','fontsize',12,'fontweight','bold');
ylabel('Relative Error','fontsize',12,'fontweight','bold');
grid on
legend('Convergence Rate','h^2');
