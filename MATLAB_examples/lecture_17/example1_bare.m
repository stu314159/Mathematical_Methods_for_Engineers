clear
clc
close 'all'

n = 150; 

f = @(x) ex1(x); 
p = pi;

a0 = (1/p)*integral(f,-p,p);
FF = @(x) a0/2;

for i = 1:n
an = (1/p)*integral(@(x) f(x).*cos(i*pi*x/p),-p,p);
bn = (1/p)*integral(@(x) f(x).*sin(i*pi*x/p),-p,p);
FF = @(x) FF(x) + an*cos(i*pi*x/p) + bn*sin(i*pi*x/p); 
end

Nx = 1000;
X = linspace(-p,p,Nx);

plot(X,f(X),'-b',...
    X,FF(X),'--r',...
    'LineWidth',3)
title_str = sprintf('Example 1, n = %d',n);
title(title_str,'FontSize',16,...
    'FontWeight','bold');
xlabel('X','FontSize',14,...
    'FontWeight','bold');
ylabel('f(X)','FontSize',14,...
    'FontWeight','bold');
grid on
legend('f(x)','FF(x)')
set(gca,'FontSize',12,...
    'FontWeight','bold');


%% Local functions 
function y = ex1(x)
[m,n] = size(x); 
assert(min(size(x))==1,'x must be a vector!')
y = nan(m,n); 
for i = 1:length(x) 
    if (x(i) >= -pi) && (x(i) < 0) 
        y(i) = 0;
    elseif (x(i) >= 0) && (x(i) <= pi)
        y(i) = pi - x(i);
    end
end
end
