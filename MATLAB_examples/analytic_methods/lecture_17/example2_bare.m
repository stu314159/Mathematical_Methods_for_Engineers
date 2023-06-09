%% example2_bare.m

clear
clc
close 'all'

n = 150; 

f = @(x) ex2(x); 
p = pi;

figure(1)
fplot(@(x) f(x), [-p,p],'-b','linewidth',3);
grid on
title('Lecture 17 Example 2','fontsize',14,'fontweight','bold');
xlabel('X','fontsize',12,'fontweight','bold');
ylabel('f(x)','fontsize',12,'fontweight','bold');
set(gca,'fontsize',12,'fontweight','bold');


a0 = (1/p)*integral(f,-p,p);
FF = @(x) a0/2;

for i = 1:n
an = (1/p)*integral(@(x) f(x).*cos(i*pi*x/p),-p,p);
fprintf('a_%d = %g \n',i,an);
bn = (1/p)*integral(@(x) f(x).*sin(i*pi*x/p),-p,p);
FF = @(x) FF(x) + an*cos(i*pi*x/p) + bn*sin(i*pi*x/p); 
end

Nx = 1000;
X = linspace(-p,p,Nx);

figure(2)
plot(X,f(X),'-b',...
    X,FF(X),'--r',...
    'LineWidth',3)
title_str = sprintf('Example 2, N = %d',n);
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
function y = ex2(x)
[m,n] = size(x); % write your function to expect vector inputs
y = nan(m,n); % construct y such that it has the same shape as x
for i = 1:length(x) %<- this implies that x will be a 1-D array
    if (x(i) > -pi) && (x(i) < 0) %<-- make sure you understand this
        y(i) = -1;
    elseif (x(i) >= 0) && (x(i) < pi)
        y(i) = 1;
    end
end
end

