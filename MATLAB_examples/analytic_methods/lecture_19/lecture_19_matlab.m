%% Lecture 19 MATLAB - Fourier-Bessel Series Expansion
% Written by: CAPT Stu Blair
% Date: 4 Oct 2020

clear
clc
close 'all'

%% Parameters
f = @(x) x; % function to be expanded
%f = @(x) x.*(x-3); % alternative function.
N = 15; % number of terms for the expansion
nu = 1; % order of Bessel function
kind = 1; % kind of Bessel function

a = 0; b = 3; % boundaries
%% Find eigenvalues

k = besselzero(nu,N,kind); % get roots
% fprintf('First %d roots of bessel function of %d kind of order %g: \n',...
%     N,nu,kind);
% disp(k);

% get the eigenvalues
alpha = k/b; % alpha*b = a root.  For this problem "3*alpha" is a root.

%% Build expansion
cn = nan(N,1);
error_norm = nan(N,1); % measure square norm of error after each term.

FB = @(x) 0;
for i = 1:N
    % compute the i-th coefficient
    cn(i) =  integral(@(x) f(x).*besselj(nu,alpha(i)*x).*x,a,b)./...
        integral(@(x) x.*besselj(nu,alpha(i)*x).^2,a,b);
    % update the Fourier-Bessel expansion
    FB = @(x) FB(x) + cn(i)*besselj(nu,alpha(i)*x);
    
    % calculate sqruare norm of the relative "error"
    err_fn = @(x) FB(x) - f(x); 
    error_norm(i) = integral(@(x) err_fn(x).^2,a,b)./...
        integral(@(x) f(x).^2,a,b); % normalize error by size of function.
end

%% Plot the result
figure(1)
fplot(FB,[a,b],'-.g','linewidth',3);
hold on
fplot(f,[a,b],'--b','linewidth',3);
hold off
grid on
xlabel('X','fontsize',14,'fontweight','bold');
ylabel('f(X)','fontsize',14,'fontweight','bold');
titlestr = sprintf('Fourier-Bessel expansion, N = %d',N);
title(titlestr,'fontsize',16,'fontweight','bold');
set(gca,'fontsize',12,'fontweight','bold');

%% Plot the error
figure(2)
f2 = semilogy(1:N,error_norm,'-ok','linewidth',3);
%f2.MarkerIndices = 1:5:length(f2.YData);
title('Convergence behavior','fontsize',16,'fontweight','bold');
grid on
xlabel('Number of Fourier-Bessel Terms','fontsize',14,'fontweight','bold');
ylabel('Relative Error','fontsize',14,'fontweight','bold');
set(gca,'fontsize',12,'fontweight','bold');