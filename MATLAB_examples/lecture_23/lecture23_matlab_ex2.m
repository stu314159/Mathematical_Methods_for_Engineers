%% Lecture 23 MATLAB Ex 2
clear
clc
close 'all'

%% Example: homogeneous Dirichlet BCs
L = 1; % length of the domain
k = 1.5; % thermal diffusivity

N = 25; % number of terms to the series solution

alpha = @(n) (2*n - 1)*pi./(2*L);
F = @(x,n) sin(alpha(n).*x);
G = @(t,n) exp(-k.*(alpha(n).^2).*t);

pick_IC = 1;
% 1 = smooth IC
% 2 = discontinuous IC
switch pick_IC
    case 1
        f = @(x) x.*(1-x);
    case 2
        f = @(x) disc_IC(x,L);
    otherwise
        error('Invalid case! \n');
end
% initialize my series solution
u = @(x,t) 0;
for n = 1:N
    % essentially doing the sine-series half-wave expansion
    % compute the coefficient
    %bn = (2/L)*integral(@(x) f(x).*F(x,n),0,L);
    %^^^ still true, but alternative implementation below...
    bn = integral(@(x) f(x).*F(x,n),0,L)./...
        integral(@(x) F(x,n).^2,0,L);
    
    % add the term to the series solution
    u = @(x,t) u(x,t) + bn*F(x,n).*G(t,n);
end
% make a discrete X-axis
Nx = 1000;
X = linspace(0,L,Nx);
%% plot the result
figure(1)
plot(X,u(X,0),'-b',...
    X,u(X,0.1),'-.g',...
    X,u(X,0.5),'--r','linewidth',3);
title('Insulated Boundary','fontsize',16,'fontweight','bold');
xlabel('X','fontsize',14,'fontweight','bold');
ylabel('u(X,t)','fontsize',14,'fontweight','bold');
grid on
set(gca,'fontsize',12,'fontweight','bold');
legend('t = 0','t = 0.1','t = 0.5');

%% time-dependent plot
Tmax = 0.5;
NT = 40;
T = linspace(0,Tmax,NT);
figure(2)
for n = 1:NT
   plot(X,u(X,T(n)),'-b','linewidth',2)
   title_str = sprintf('Lecture #23 Example, t = %g ',T(n));
   title(title_str,'fontsize',16,'fontweight','bold');
   xlabel('X','fontsize',14,'fontweight','bold');
   ylabel('u(X,T)','fontsize',14,'fontweight','bold');
   grid on
   set(gca,'fontsize',12,'fontweight','bold');
   axis([0 L -0.2 1.2]);
   pause(Tmax/(NT-1)*10);
end

%% Local functions
function y = disc_IC(x,L)
[m,n] = size(x);
y = nan(m,n);
for i = 1:length(x)
   if (x(i) > 0) && (x(i) < L/4)
       y(i) = x(i);
   elseif(x(i) >= L/4) && (x(i) < L/2)
       y(i) = 1;
   elseif(x(i) >= L/2) && (x(i) < 3*L/4)
       y(i) = 0;
   elseif(x(i) >= 3*L/4) && (x(i) < L)
       y(i) = L - x(i);
   end
end
end