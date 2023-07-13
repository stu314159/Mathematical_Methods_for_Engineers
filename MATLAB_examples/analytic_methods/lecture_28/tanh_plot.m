%%tanh_plot.m

clear
clc
close 'all'

f = @(x) tanh(x);

a = 0; b = 10;

figure(1)
fplot(f,[a b],'linewidth',3)
title('Plot of $$\tanh{(\nu)}$$','interpreter','latex');
xlabel('\nu');
ylabel('$$\tanh{(\nu)}$$','Interpreter','latex');
grid on
set(gca,'fontsize',12,'fontweight','bold');
