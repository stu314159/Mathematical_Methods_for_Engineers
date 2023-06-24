%% Lecture 27 MATLAB example
clear
clc
close 'all'

%% Define parameters
L = 1;
k = 0.5; % thermal diffusivity
r = 150; % constant heat production in the domain
u1 = 100; % right-hand boundary condition (time independent)
f = @(x)  300*x.*(1-x); % initial condition
alpha = @(n) n*pi;
F = @(x,n) sin(alpha(n).*x);
G = @(t,n) exp(-k*(alpha(n).^2).*t);

%% Construct truncated series solution
N = 100; % number of terms to use
V = @(x,t) 0;
psi = @(x) -r./(2*k).*(x.^2) + (u1 + r/(2*k)).*x;

for n = 1:N
   % compute the nth Fourier Coefficient
   bn = 2/L*integral(@(x) (f(x)-psi(x)).*F(x,n),0,L);
   % add the next term to the series
   V = @(x,t) V(x,t) + bn*F(x,n).*G(t,n);
end

U = @(x,t) V(x,t)+psi(x);

%% Plot the solution at fixed times
Nx = 1000;
X = linspace(0,L,Nx);
figure(1)
plot(X,U(X,0),'-b',...
    X,U(X,0.1),'--g',...
    X,U(X,inf),':k','linewidth',3);
%title('Lecture 27 Example u(x,t)','fontsize',16,'fontweight','bold');
title('Lecture 27 Example u(x,t)');

xlabel('X','fontsize',14,'fontweight','bold');
ylabel('U(X,T)','fontsize',14,'fontweight','bold');
grid on
set(gca,'fontsize',12,'fontweight','bold');
legend('t=0','t=0.1','t=\infty','location','best');
%axis([0 1 0 4])
% %% Time-dependent plot
% NT = 50;
% Tmax = 2;
% T = linspace(0,Tmax,NT);
% figure(2)
% for t = 1:NT
%     plot(X,U(X,T(t)),'linewidth',3);
%     title_str = sprintf('Lecture 27 Example: t = %g ',T(t));
%     title(title_str,'fontsize',16,'fontweight','bold');
%     xlabel('X','fontsize',14,'fontweight','bold');
%     ylabel('U(X,T)','fontsize',14,'fontweight','bold');
%     grid on
%     set(gca,'fontsize',12,'fontweight','bold');
%     axis([0 1 -2 10]);
%     pause(Tmax/(NT-1));    
% end
% hold on
% plot(X,psi(X),'--r','linewidth',3);
% hold off
% legend('U(X,\infty)','\psi(X)','location','best');

%% Fixed in time plot
figure(3)
plot(X,psi(X),'-b','linewidth',3)
title('Lecture 27 Example: $$\psi(x)$$',...
    'Interpreter','latex')
grid on
xlabel('X','fontsize',14,'fontweight','bold');
ylabel('\psi(X)','fontsize',14,'fontweight','bold');
set(gca,'fontsize',12,'fontweight','bold');
%axis([0 1 0 4])