clear
clc 
close 'all'

f_even = @(x) cosh(x);

p = 5;
xMin = -p; xMax = p; nX = 1000;
X = linspace(xMin,xMax,nX);

figure(1)
plot(X,f_even(X),'-b','linewidth',3);
title('Even Function','fontsize',16,...
    'FontWeight','bold');
grid on;
xlabel('X','fontsize',14,'FontWeight','bold');
ylabel('f(X)','FontSize',14,'FontWeight','bold');
set(gca,'fontsize',12,'fontweight','bold');

axh = gca; % use current axes
color = 'k'; % black, or [0 0 0]
linestyle = '--'; % dotted
line(get(axh,'XLim'), [0 0], ...
    'Color', color, ...
    'LineStyle', linestyle,...
    'LineWidth',3);
line([0 0], get(axh,'YLim'), ...
    'Color', color, ...
    'LineStyle', linestyle,...
    'LineWidth',3);


f_odd = @(x) x.^3;
figure(2)
plot(X,f_odd(X),'-b','linewidth',3);
title('Odd Function','fontsize',16,...
    'FontWeight','bold');
grid on;
xlabel('X','fontsize',14,'FontWeight','bold');
ylabel('f(X)','FontSize',14,'FontWeight','bold');
set(gca,'fontsize',12,'fontweight','bold');
axh = gca; % use current axes
line(get(axh,'XLim'), [0 0], ...
    'Color', color, ...
    'LineStyle', linestyle,...
    'LineWidth',3);
line([0 0], get(axh,'YLim'), ...
    'Color', color, ...
    'LineStyle', linestyle,...
    'LineWidth',3);