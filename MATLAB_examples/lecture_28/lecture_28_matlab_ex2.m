%% Lecture 28 MATLAB example 2
clear
clc
close 'all'

%% Parameters
alpha_sq = 5;
L = 1;

N = 30;

nu = @(n) (2*n-1)*pi/2;
Fn = @(x,n) sin(nu(n).*x);
Gn = @(t,n) cos(alpha_sq*nu(n).*t);
f = @(x) x; % initial displacement
u = @(x,t) 0; % initial velocity

for n = 1:N
   ef_mag = integral(@(x) Fn(x,n).^2,0,L);
   an = integral(@(x) f(x).*Fn(x,n),0,L)./ef_mag;
   
   u = @(x,t) u(x,t) + an.*Fn(x,n).*Gn(t,n);
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
    title_str = sprintf('Lecture 28 Example 2, t = %g ',T(t));
    title(title_str,'fontsize',16,'fontweight','bold');
    xlabel('X','fontsize',14,'fontweight','bold');
    ylabel('U(X,T) Angular Displacement', 'fontsize',14,...
        'fontweight','bold');
    grid on
    set(gca,'fontsize',12,'fontweight','bold');
    axis([0 L -L L]);
    pause(Tmax/(Nt-1));    
end