clear
clc
close 'all'

N = 500; % number of eigenvalues
a = 0; b = 3; % bounds of the domain
nu = 1; kind = 1;
k = besselzero(nu,N,kind); % get roots
alpha = k/b;

f = @(x) x; 
cn = nan(N,1); % store the coefficients (optional)
rel_err = nan(N,1); 

FB = @(x) 0; % initialize the Fourier-Bessel expansion
for n = 1:N
    % compute the i-th coefficient
    cn(n) =...
        integral(@(x) f(x).*besselj(nu,alpha(n)*x).*x,a,b)./...
        integral(@(x) x.*besselj(nu,alpha(n)*x).^2,a,b);
    % update the Fourier-Bessel expansion
    FB = @(x) FB(x) + cn(n)*besselj(nu,alpha(n)*x);

    % calculate sqruare norm of the relative "error"
    err_fn = @(x) FB(x) - f(x); 
    rel_err(n) = integral(@(x) err_fn(x).^2,a,b)./...
        integral(@(x) f(x).^2,a,b); 
end

figure(1)
loglog(1:N,rel_err,'-b',...
    'LineWidth',3);
title('\textbf{Convergence of Fourier-Bessel Expansion of } $$f(x)=x$$',...
    'Interpreter','latex');
ylabel('Relative Error','FontSize',14,...
    'FontWeight','bold');
xlabel('Number of Terms','FontSize',14,...
    'FontWeight','bold');
grid on
set(gca,'FontSize',12,'FontWeight','bold');
