clear
clc
close 'all'

f = @(x) sin(x) + x.*cos(x);

xMin = 0; xMax = 25;
Nx = 1000;
X = linspace(xMin,xMax,Nx);

figure(1)
plot(X,f(X),'-b','linewidth',3)
title('$$\sin{\alpha} + \alpha \cos{\alpha}$$',...
    'interpreter','latex',...
    'FontSize',18,'FontWeight','bold');
grid on
xlabel('\alpha','fontsize',16,'fontweight','bold');
set(gca,'fontsize',12,'fontweight','bold');

axh = gca; % use current axes
color = 'k'; % black, or [0 0 0]
linestyle = '--'; % dotted
line(get(axh,'XLim'), [0 0], ...
    'Color', color, ...
    'LineStyle', linestyle,...
    'LineWidth',3);

%%
f = @(x) x + tan(x);

figure(2)
plot(X,f(X),'-b','linewidth',3)
title('$$\alpha + \tan{\alpha}$$',...
    'interpreter','latex',...
    'FontSize',18,'FontWeight','bold');
grid on
xlabel('\alpha','fontsize',16,'fontweight','bold');
set(gca,'fontsize',12,'fontweight','bold');
axis([xMin xMax -40 40])
axh = gca; % use current axes
color = 'k'; % black, or [0 0 0]
linestyle = '--'; % dotted
line(get(axh,'XLim'), [0 0], ...
    'Color', color, ...
    'LineStyle', linestyle,...
    'LineWidth',3);

figure (3)
plot(X,tan(X),'-b',...
    X,-X,'--r','linewidth',3)
% title('$$\sin{\alpha} + \alpha \cos{\alpha}$$',...
%     'interpreter','latex',...
%     'FontSize',18,'FontWeight','bold');
grid on
xlabel('\alpha','fontsize',16,'fontweight','bold');
set(gca,'fontsize',12,'fontweight','bold');
axis([xMin xMax -40 40])
axh = gca; % use current axes
color = 'k'; % black, or [0 0 0]
linestyle = '--'; % dotted
line(get(axh,'XLim'), [0 0], ...
    'Color', color, ...
    'LineStyle', linestyle,...
    'LineWidth',3);
legend('tan(\alpha)','-\alpha')