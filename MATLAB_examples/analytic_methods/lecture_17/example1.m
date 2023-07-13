%% Example 1
clear
clc 
close 'all'

%% Example 1
f = @(x) ex1(x); %<- not mandatory but provides consistent syntax

N1 = 5;
p = pi;

a0 = (1/p)*integral(f,-p,p);
FF1 = @(x) a0/2;

for i = 1:N1
    an = (1/p)*integral(@(x) f(x).*cos(i*pi*x/p),-p,p);
    bn = (1/p)*integral(@(x) f(x).*sin(i*pi*x/p),-p,p);
    
    FF1 = @(x) FF1(x) + an*cos(i*pi*x/p) + bn*sin(i*pi*x/p);
end
FF2 = @(x) FF1(x);
N2 = 150;
for i = (N1+1):N2
    an = (1/p)*integral(@(x) f(x).*cos(i*pi*x/p),-p,p);
    bn = (1/p)*integral(@(x) f(x).*sin(i*pi*x/p),-p,p);
    
    FF2 = @(x) FF2(x) + an*cos(i*pi*x/p) + bn*sin(i*pi*x/p);
end

figure(1)
fplot(@(x) f(x), [-p,p],'-b','linewidth',3);
grid on
title('Lecture 17 Example 1','fontsize',14,'fontweight','bold');
xlabel('X','fontsize',12,'fontweight','bold');
ylabel('f(x)','fontsize',12,'fontweight','bold');
set(gca,'fontsize',12,'fontweight','bold');

figure(2)
fplot(@(x)f(x),[-p,p],'-b','linewidth',2)
hold on
fplot(@(x) FF1(x),[-p,p],'-.g','linewidth',2)
fplot(@(x) FF2(x),[-p,p],'-or','linewidth',2)
hold off
grid on
title('Lecture 17 Example 1','fontsize',14,'fontweight','bold');
xlabel('X','fontsize',12,'fontweight','bold');
ylabel('F(x)','fontsize',12,'fontweight','bold');
set(gca,'fontsize',12,'fontweight','bold');
legend('f(x)','N=5','N=15','location','best');


%% Local functions for examples
function y = ex1(x)
[m,n] = size(x); % write your function to expect vector inputs
y = nan(m,n); % construct y such that it has the same shape as x
for i = 1:length(x) %<- this implies that x will be a 1-D array
    if (x(i) > -pi) && (x(i) < 0) %<-- make sure you understand this
        y(i) = 0;
    elseif (x(i) >= 0) && (x(i) < pi)
        y(i) = pi - x(i);
    end
end
end
