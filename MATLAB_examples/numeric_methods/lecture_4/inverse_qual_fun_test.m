% inverse quadratic function.
clear
clc
close 'all'

f = @(x) x.^2 - 2;

a = 0; b = 1.45; c = 2;


% x(y) is the inverse quadratic interpolant of f(x)
x = @(y) (y - f(a)).*(y-f(b)).*c./((f(c)-f(a)).*(f(c)-f(b))) + ...
    (y - f(b)).*(y - f(c)).*a./((f(a)-f(b)).*(f(a)-f(c))) + ...
    (y-f(c)).*(y-f(a)).*b./((f(b)-f(c)).*(f(b)-f(a)));
% root is approximately at x(0)
% it doesn't matter how a,b, and c are ordered; only
% that the root lies somewhere in the real line spanned by
% those three points.

X = linspace(0,2,500);
Y = linspace(-2,2,500);

figure (1)
plot(X,f(X),'-.g',...
    x(Y),Y,'-r','linewidth',2)
title('Inverse Quadratic Interpolation','fontsize',14,...
    'fontweight','bold');
ylabel('Y','fontsize',14,'fontweight','bold');
xlabel('X','fontsize',14,'fontweight','bold');
legend('f(x)','Interpolant','location','best');
grid on
set(gca,'fontsize',12,'fontweight','bold');
