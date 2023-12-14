%gauss_quad_convergence.m

clear
clc
close all

f = @(x) exp(-x.^2);
a = -3;
b = 3;

intF_exact = integral(f,a,b);

N = 2:20;

relError = NaN(1,length(N));

for p = N
   intF = myGaussQuad1D(f,a,b,p);
   relError(p-1) = norm(intF-intF_exact,2)/intF_exact;
end

semilogy(N,relError,'linewidth',2);
title('Convergence of Gauss Quadrature','fontsize',16,'fontweight','bold');
xlabel('Number of Gauss Points (P)','fontsize',14,'fontweight','bold');
ylabel('Relative Error','fontsize',14,'fontweight','bold');
set(gca,'fontsize',12,'fontweight','bold');
grid on

loglog(N,relError,'linewidth',2);
title('Convergence of Gauss Quadrature','fontsize',16,'fontweight','bold');
xlabel('Number of Gauss Points (P)','fontsize',14,'fontweight','bold');
ylabel('Relative Error','fontsize',14,'fontweight','bold');
set(gca,'fontsize',12,'fontweight','bold');
grid on