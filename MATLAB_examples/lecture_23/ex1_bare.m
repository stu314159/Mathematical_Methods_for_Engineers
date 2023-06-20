clear
clc
close 'all'

%% Example: homogeneous Dirichlet BCs
L = 1; % length of the domain
alpha_sq = 0.1; % thermal diffusivity

N = 100; % number of terms to the series solution

F = @(x,n) sin(n.*pi.*x./L);
G = @(t,n) exp(-((n.*pi./L).^2)*alpha_sq.*t);


pick_IC = 2;
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

%% Plot the initial condition
% make a discrete X-axis
Nx = 1000;
X = linspace(0,L,Nx);

figure(1)
plot(X,f(X),'-c','LineWidth',3)
title('Example Initial Condition',...
    'FontSize',14,'fontweight','bold');
xlabel('X','FontSize',12,'FontWeight','bold');
set(gca,'fontsize',10,'fontweight','bold');
grid on


%% Compute the solution
% initialize my series solution
u = @(x,t) 0;
for n = 1:N
    % essentially doing the sine-series half-wave expansion
    % compute the coefficient
    cn = (2/L)*integral(@(x) f(x).*F(x,n),0,L);
    
    % add the term to the series solution
    u = @(x,t) u(x,t) + cn*F(x,n).*G(t,n);
end

%% plot the result
figure(1)
plot(X,u(X,0),'-b',...
    X,u(X,0.1),'-.g',...
    X,u(X,0.5),'--r','linewidth',3);
title_str = sprintf('Heat Equation Example, N=%d',N);
title(title_str,'FontSize',16,'FontWeight','bold');
%title('Lecture #23 Example','fontsize',16,'fontweight','bold');
xlabel('X','fontsize',14,'fontweight','bold');
ylabel('u(X,t)','fontsize',14,'fontweight','bold');
grid on
set(gca,'fontsize',12,'fontweight','bold');
legend('t = 0','t = 0.1','t = 0.5');

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