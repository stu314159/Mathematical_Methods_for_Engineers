%% Lecture 20 Matlab Example
% Written by: CAPT Stu Blair
% Date: 4 Oct 2020
clear
clc
close 'all'

%% Parameters
%f = @(x) ex1(x);
f = @(x) cos(x);

N = 12; % number of terms 
a = -1; b = 1; % boundaries

c0 = (1/2)*integral(@(x) f(x),a,b);
cn = nan(N-1,1);
error_norm = nan(N,1); % measure square norm of error after each term.

FL = @(x) c0;
err_fn = @(x) FL(x) - f(x);
error_norm(1) = integral(@(x) err_fn(x).^2,a,b)./...
        integral(@(x) f(x).^2,a,b); % normalize error by size of function.

for i = 1:(N-1)
    % compute the i'th coefficient
    cn(i) = ((2*i+1)/2)*integral(@(x) f(x).*legendreP(i,x),a,b);
    FL = @(x) FL(x) + cn(i)*legendreP(i,x); %update the expansion
    
    % compute the error.
    err_fn = @(x) FL(x) - f(x);
    error_norm(i+1) = integral(@(x) err_fn(x).^2,a,b)./...
        integral(@(x) f(x).^2,a,b); % normalize error by size of function.
end

%% Plot the result
figure(1)
fplot(FL,[a,b],'-.g','linewidth',3);
hold on
fplot(f,[a,b],'--b','linewidth',3);
hold off
grid on
xlabel('X','fontsize',14,'fontweight','bold');
ylabel('f(X)','fontsize',14,'fontweight','bold');
titlestr = sprintf('Fourier-Legendre expansion, N = %d',N);
title(titlestr,'fontsize',16,'fontweight','bold');
set(gca,'fontsize',12,'fontweight','bold');

%% Plot the error
figure(2)
semilogy(1:N,error_norm,'-ok','linewidth',3);
title('Convergence behavior','fontsize',16,'fontweight','bold');
grid on
xlabel('Number of Fourier-Legendre Terms','fontsize',14,'fontweight','bold');
ylabel('Relative Error','fontsize',14,'fontweight','bold');
set(gca,'fontsize',12,'fontweight','bold');
%% Local function to implement f(x)
function y = ex1(x)
[m,n] = size(x); % write your function to expect vector inputs
assert(min(m,n) == 1,'Bad input for ex1');
y = nan(m,n); % construct y such that it has the same shape as x
for i = 1:length(x) %<- this implies that x will be a 1-D array
    if (x(i) > -1) && (x(i) < 0) %<-- make sure you understand this
        y(i) = 0;
    elseif (x(i) >= 0) && (x(i) < 1)
        y(i) = 1;
    end
end
end