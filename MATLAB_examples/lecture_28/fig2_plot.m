%% fig2_plot

clear
clc
close 'all'

f1 = @(x) tan(x);

h = 1;
f2 = @(x) -x/h;

a = 1e-8; b = 15;
Nx = 1000;
X = linspace(a,b,Nx);
figure(1)
plot(X,f1(X),'-b',...
    X,f2(X),'--k','linewidth',3)
grid on
ylabel('f(x)','fontsize',14,'fontweight','bold');
xlabel('\nu','FontSize',14,'FontWeight','bold');
legend('$$\tan(\nu)$$','$$-\frac{\nu}{h}$$','interpreter','latex');

axis([a b -20 20])