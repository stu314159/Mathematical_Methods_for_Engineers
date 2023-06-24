%% Lecture 28 MATLAB
clear
clc
close 'all'

%% Parameters
k = 0.01; % thermal diffusivity
h = 0.5; % convective heat transfer coefficient

L = 1;
f = @(x) 0.*x + 1; % initial condition

N = 25; % number of eigenfunctions

ef_fun = @(x) tan(x) + x./h;


%% plot the condition for the eigenfunction
figure(1)
fplot(ef_fun,[0,20],'-b','linewidth',3);
xlabel('alpha','fontsize',14,'fontweight','bold');
ylabel('f(alpha)','fontsize',14,'fontweight','bold');
axis([0 20 -20 20]);
grid on
set(gca,'fontsize',12,'fontweight','bold');

%% Find the desired number of eigenvalues
alpha = nan(N,1);
x0 = [pi/2+1e-5,3*pi/2-1e-5];
for i = 1:N
   [x,fval,exit_flag] = fzero(ef_fun,x0);
   alpha(i) = x;
   x0 = x0 + pi;
end
% do some crude error checking
assert(min(diff(alpha))>0,...
    'Error! Something is wrong with your eigenvalues!');

%% Construct solution
Fn = @(x,n) sin(alpha(n).*x);
Gn = @(t,n) exp(-(alpha(n).^2).*k.*t); %<--bug fixed!

u = @(x,t) 0;

for n = 1:N
   ef_mag = integral(@(x) Fn(x,n).*Fn(x,n),0,L);
   Bn = integral(@(x) f(x).*Fn(x,n),0,L)./ef_mag;
   
   u = @(x,t) u(x,t) + Bn*Fn(x,n).*Gn(t,n);
end

%% Plot the solution
Nx = 1000;
X = linspace(0,L,Nx);

Tmax = 5;
Nt = 50;
T = linspace(0,Tmax,Nt);

figure(2)
for t = 1:Nt
    plot(X,u(X,T(t)),'-b','linewidth',3);
    title_str = sprintf('Lecture 28 Example, t = %g ',T(t));
    title(title_str,'fontsize',16,'fontweight','bold');
    grid on
    set(gca,'fontsize',12,'fontweight','bold');
    axis([0 L -0.05 1.5]);
    pause(Tmax/(Nt-1));
    
end
