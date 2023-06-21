%% Lecture 24 MATLAB
clear
clc
close 'all'

%% Example Problem
L = 3;
alpha_sq = 1;% T/rho
alpha = sqrt(alpha_sq);
N = 50;

ex_select = 2;
switch ex_select
    case 1
        f = @(x) ex1(x,L);
    case 2
        f = @(x) ex2(x,L);
    case 3
        alpha = 1;
        f = @(x) exp(-20*(x-1).^2) + exp(-20*(x-2).^2);
    case 4
        alpha = 1;
        f = @(x) exp(-20*(x-0.5).^2) + exp(-20*(x-1.5).^2) + ...
            exp(-20*(x-2.5).^2);
    otherwise
        error('Invalid case!');
end

g = @(x) x.*0;
%g = @(x) -x.^2;


u = @(x,t) 0;

for n = 1:N
    % compute an
    an = (2/L)*integral(@(x) f(x).*sin(n*pi*x./L),0,L);
    % compute bn
    bn = ...
        (2/(alpha*n*pi))*...
        integral(@(x) g(x).*sin(n*pi*x./L),0,L);
    
    % update the approximate solution
    u = @(x,t) u(x,t) + ...
        (an*cos(alpha.*n*pi*t./L) + ...
        bn*sin(alpha.*n*pi*t./L)).*sin(n*pi*x./L); 
end
%% make discrete space and time spaces
Tmax = 3;
NT = 50;
T = linspace(0,Tmax,NT);

Nx = 500;
X = linspace(0,L,Nx);

%% time-dependent plot

figure(1)
for n = 1:NT
   plot(X,u(X,T(n)),'-b','linewidth',2); 
   title_str = sprintf('Lecture 24 Example, t = %g ',T(n));
   title(title_str,'fontsize',16,'fontweight','bold');
   xlabel('X','fontsize',14,'fontweight','bold');
   ylabel('u(X,T)','fontsize',14,'fontweight','bold');
   grid on
   set(gca,'fontsize',12,'fontweight','bold');
   axis([0 L -1.5 1.5]);
   pause(Tmax/(NT-1));
end

%% fixed plot, multiple data series
figure(2)
plot(X,u(X,0),'-b',...
    X,u(X,1.0),'-.g',...
    X,u(X,2.0),'--r',...
    X,u(X,3.0),':k','linewidth',3);
title('Lecture 25 Wave Equation Example',...
    'fontsize',16,'fontweight','bold');
xlabel('X','fontsize',14,'fontweight','bold');
ylabel('u(X,T)','fontsize',14,...
    'fontweight','bold');
grid on
set(gca,'fontsize',12,'fontweight','bold');
axis([0 L -1.0 1.0]);
legend('t = 0','t = 1.0','t = 2.0','t = 3.0',...
    'location','best');

%% Fixed Plot, single time step
ts = 3.0;
figure(3)
plot(X,u(X,ts),'-b','Linewidth',3);
title_str = ...
    sprintf('Lecture 24 Example, t = %g',ts);
title(title_str,'fontsize',16,'fontweight','bold');
xlabel('X','fontsize',14,'fontweight','bold');
ylabel('u(X,T)','fontsize',14,'fontweight','bold');
grid on
set(gca,'fontsize',12,'fontweight','bold');
axis([0 L -1 1]);


%% Local functions
function y = ex1(x,L)
[m,n] = size(x);
y = nan(m,n);
for i = 1:length(x)
   if (x(i)>0)&& (x(i) < L/2)
       y(i) = (2/3).*x(i);
   elseif(x(i) >= L/2) && (x(i)<L)
       y(i) = (2/3)*(L - x(i));
   end
end
end

function y = ex2(x,L)
[m,n] = size(x);
y = nan(m,n);
for i = 1:length(x)
    if(x(i)>0) && (x(i)< L/3)
        y(i) = 0;
    elseif(x(i) >= L/3) && (x(i)< L/2)
        y(i) = (x(i)-L/3);
    elseif(x(i) >= L/2)&&(x(i) < (2*L/3))
        y(i) = L/6 - (x(i)-L/2);
    elseif(x(i) >= 2*L/3)&&(x(i) < L)
        y(i) = 0;
    end
end
end