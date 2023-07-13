%% Lecture 35 example
clear
clc
close 'all'

%% Parameters
N = 50;

c = 1; % radius of sphere
Psi = @(r) 100.*r;
R = @(r,n) sin(n.*pi.*r);
T = @(t,n) exp(-((n.*pi).^2).*t);

V = @(r,t) 0; % initialize my approximation for V

U_gold = @(r,t) 100;

for n = 1:N
   % compute coefficient
   cn = 2*integral(@(r) -Psi(r).*R(n,r),0,c);
   
   % update solution to V
   V = @(r,t) V(r,t) + cn.*R(r,n).*T(t,n);
   
   % update given solution (from textbook...)
   U_gold = @(r,t) U_gold(r,t) + ...
       (200./(pi.*r)).*((-1).^n)./(n).*R(r,n).*T(t,n);
end

% construct solution U
U = @(r,t) (1./r).*(V(r,t)+Psi(r));

%% Time-dependent Plot
Tmax = 1;
NT = 50;
T = linspace(0,Tmax,NT);

NR = 250;
R = linspace(0,c,NR);

% figure(1)
% for i = 1:NT
%    plot(R,U(R,T(i)),'-k',...
%        R,U_gold(R,T(i)),'--r','linewidth',3);
%    title_str = sprintf('Lecture 35 Example, t = %g',T(i));
%    title(title_str,'fontsize',16,'fontweight','bold');
%    xlabel('R','fontsize',14,'fontweight','bold');
%    ylabel('U(R,T)','fontsize',14,'fontweight','bold');
%    grid on
%    set(gca,'fontsize',12,'fontweight','bold');
%    legend('U','U_{gold}','location','southeast');
%    axis([0 c -5 105]);
%    pause(5*Tmax/(NT-1));    
% end

%% Plots at discrete times

t = 0.3;
figure(2)
plot(R,U(R,t),'-k','linewidth',3)
grid on
title_str = sprintf('Lecture 35 example, t = %g',t);
title(title_str,'fontsize',16,'fontweight','bold');
xlabel('R','fontsize',14,'fontweight','bold');
ylabel('U(R,T)','fontsize',14,'fontweight','bold');
set(gca,'fontsize',12,'fontweight','bold');
axis([0 c -5 105]);
