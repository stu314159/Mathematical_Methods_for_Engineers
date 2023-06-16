clear
clc
close 'all'

f = @(x) ex1(x);

N = 15; % number of terms 
a = -1; b = 1; % boundaries

% handle Po coefficient separately
c0 = (1/2)*integral(@(x) f(x),a,b);
cn = nan(N-1,1);
error_norm = nan(N,1); 

FL = @(x) c0;

% calculate relative error
err_fn = @(x) FL(x) - f(x);
error_norm(1) = integral(@(x) err_fn(x).^2,a,b)./...
        integral(@(x) f(x).^2,a,b); 

for n = 1:(N-1)
    % compute the n'th coefficient
    cn(n) = ((2*n+1)/2)*integral(@(x) f(x).*legendreP(n,x),a,b);
    FL = @(x) FL(x) + cn(n)*legendreP(n,x); %update the expansion
    
    % compute the error.
    err_fn = @(x) FL(x) - f(x);
    error_norm(n+1) = integral(@(x) err_fn(x).^2,a,b)./...
        integral(@(x) f(x).^2,a,b); % normalize error by size of function.
end

%% Plot the result
Nx = 1000;
X = linspace(a,b,Nx);

figure(1)
plot(X,FL(X),'-g',...
    X,f(X),'--b',...
    'LineWidth',3);
grid on
xlabel('X','fontsize',14,'fontweight','bold');
ylabel('f(X)','fontsize',14,'fontweight','bold');
titlestr = ...
    sprintf('Fourier-Legendre expansion, N = %d',N);
title(titlestr,'fontsize',16,'fontweight','bold');
set(gca,'fontsize',12,'fontweight','bold');

%% Plot the error
figure(2)
loglog(1:N,error_norm,'-ok','linewidth',3);
title('Convergence behavior','fontsize',16,'fontweight','bold');
grid on
xlabel('Number of Fourier-Legendre Terms','fontsize',14,'fontweight','bold');
ylabel('Relative Error','fontsize',14,'fontweight','bold');
set(gca,'fontsize',12,'fontweight','bold');


%% Local functions
function y = ex1(x)
[m,n] = size(x); 
% expect vector inputs.
assert(min(m,n) == 1,'Bad input for ex1');
% construct y so that it has the same shape as x
y = nan(m,n); 
for i = 1:length(x) 
    if (x(i) > -1) && (x(i) < 0) 
        y(i) = 0;
    elseif (x(i) >= 0) && (x(i) < 1)
        y(i) = 1;
    end
end
end