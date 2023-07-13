clear
clc
close 'all'

N = 150; % number of eigenvalues
a = 0; b = 3; % bounds of the domain
nu = 1; kind = 1;
k = besselzero(nu,N,kind); % get roots
alpha = k/b;

f = @(x) x; 
cn = nan(N,1); % store the coefficients (optional)

FB = @(x) 0; % initialize the Fourier-Bessel expansion
for n = 1:N
    % compute the i-th coefficient
    cn(n) =...
        integral(@(x) f(x).*besselj(nu,alpha(n)*x).*x,a,b)./...
        integral(@(x) x.*besselj(nu,alpha(n)*x).^2,a,b);
    % update the Fourier-Bessel expansion
    FB = @(x) FB(x) + cn(n)*besselj(nu,alpha(n)*x);
end

Nx = 1000;
X = linspace(a,b,Nx);

figure(1)
plot(X,FB(X),'-b','LineWidth',3);
xlabel('X','fontsize',14,'fontweight','bold');
ylabel('f(X)','fontsize',14,'fontweight','bold');
titlestr = sprintf('Fourier-Bessel expansion, N = %d',N);
title(titlestr,'fontsize',16,'fontweight','bold');
grid on
set(gca,'fontsize',12,'fontweight','bold');
