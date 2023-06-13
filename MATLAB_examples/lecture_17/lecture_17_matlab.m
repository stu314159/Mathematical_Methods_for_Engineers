%% Lecture 17 Fourier Series in MATLAB
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

%% Example 2
clear;

f = @(x) ex2(x); %<- not mandatory but provides consistent syntax

N1 = 5;
p = pi;
FF1 = @(x) 0;

for i = 1:N1
    bn = (2/p)*integral(@(x) f(x).*sin(i*pi*x/p),0,p); %notice bounds of integration.
    
    FF1 = @(x) FF1(x) + bn*sin(i*pi*x/p);
end
FF2 = @(x) FF1(x);
N2 = 15;
for i = (N1+1):N2    
    bn = (2/p)*integral(@(x) f(x).*sin(i*pi*x/p),0,p);
    
    FF2 = @(x) FF2(x) + bn*sin(i*pi*x/p);
end

figure(2)
fplot(@(x)f(x),[-p,p],'-b','linewidth',2)
hold on
fplot(@(x) FF1(x),[-p,p],'-.g','linewidth',2)
fplot(@(x) FF2(x),[-p,p],'-or','linewidth',2)
hold off
grid on
title('Lecture 17 Example 2 (Sine Series)','fontsize',14,'fontweight','bold');
xlabel('X','fontsize',12,'fontweight','bold');
ylabel('F(x)','fontsize',12,'fontweight','bold');
set(gca,'fontsize',12,'fontweight','bold');
legend('f(x)','N=5','N=15','location','best');

%% Example 2 as a full Fourier Series
clear;

f = @(x) ex2(x); %<- not mandatory but provides consistent syntax

N1 = 5;
p = pi;

a0 = (1/p)*integral(@(x) f(x),-p,p);
FF1 = @(x) a0/2;

% save coefficient values for a_n and b_n
a_n = nan(16,1);
b_n = nan(16,1);
a_n(1) = a0;

for i = 1:N1
    an = (1/p)*integral(@(x) f(x).*cos(i*pi*x/p),-p,p);
    bn = (1/p)*integral(@(x) f(x).*sin(i*pi*x/p),-p,p);
    
    FF1 = @(x) FF1(x) + an*cos(i*pi*x/p) + bn*sin(i*pi*x/p);
    
    a_n(i+1) = an;
    b_n(i+1) = bn;
end
FF2 = @(x) FF1(x);
N2 = 15;
for i = (N1+1):N2
    an = (1/p)*integral(@(x) f(x).*cos(i*pi*x/p),-p,p);
    bn = (1/p)*integral(@(x) f(x).*sin(i*pi*x/p),-p,p);
    
    FF2 = @(x) FF2(x) + an*cos(i*pi*x/p) + bn*sin(i*pi*x/p);
    
    a_n(i+1) = an;
    b_n(i+1) = bn;
end

figure(3)
fplot(@(x)f(x),[-p,p],'-b','linewidth',2)
hold on
fplot(@(x) FF1(x),[-p,p],'-.g','linewidth',2)
fplot(@(x) FF2(x),[-p,p],'-or','linewidth',2)
hold off
grid on
title('Example 2 Full Fourier Series','fontsize',14,'fontweight','bold');
xlabel('X','fontsize',12,'fontweight','bold');
ylabel('F(x)','fontsize',12,'fontweight','bold');
set(gca,'fontsize',12,'fontweight','bold');
legend('f(x)','N=5','N=15','location','best');
fprintf('coefficient values for a_n: \n');
disp(a_n);
fprintf('coefficient values for b_n: \n');
disp(b_n);


%% Example 2 as a Fourier Sine Series
clear;

f = @(x) ex2(x); %<- not mandatory but provides consistent syntax

N1 = 5;
p = pi;



% save coefficient values for a_n and b_n
b_n = nan(16,1);
FF1 = @(x) 0;

for i = 1:N1
    bn = (2/p)*integral(@(x) f(x).*sin(i*pi*x/p),0,p);
    
    FF1 = @(x) FF1(x) + bn*sin(i*pi*x/p);
    
    b_n(i+1) = bn;
end
FF2 = @(x) FF1(x);
N2 = 15;
for i = (N1+1):N2
   % an = (1/p)*integral(@(x) f(x).*cos(i*pi*x/p),-p,p);
    bn = (1/p)*integral(@(x) f(x).*sin(i*pi*x/p),-p,p);
    
    FF2 = @(x) FF2(x) + bn*sin(i*pi*x/p);
    
    %a_n(i+1) = an;
    b_n(i+1) = bn;
end

figure(3)
fplot(@(x)f(x),[-p,p],'-b','linewidth',2)
hold on
fplot(@(x) FF1(x),[-p,p],'-.g','linewidth',2)
fplot(@(x) FF2(x),[-p,p],'-or','linewidth',2)
hold off
grid on
title('Example 2 Fourier Sine Series','fontsize',14,'fontweight','bold');
xlabel('X','fontsize',12,'fontweight','bold');
ylabel('F(x)','fontsize',12,'fontweight','bold');
set(gca,'fontsize',12,'fontweight','bold');
legend('f(x)','N=5','N=15','location','best');
fprintf('coefficient values for a_n: \n');
%disp(a_n);
fprintf('coefficient values for b_n: \n');
disp(b_n);


%% Example 3 part 1: cosine series expansion
clear;

f = @(x) x.^2;
L = 2; %<--to plot we need to define this
p = L;
N1 = 3; 
a0 = (2/p)*integral(@(x)f(x),0,p);

FF1 = @(x) a0/2;

for i = 1:N1
    an = (2/p)*integral(@(x) f(x).*cos(i*pi*x/p),0,p);
    
    FF1 = @(x) FF1(x) + an*cos(i*pi*x/p); 
end
FF2 = @(x) FF1(x);
N2 = 5;
for i = (N1+1):N2
    an = (2/p)*integral(@(x) f(x).*cos(i*pi*x/p),0,p);
        
    FF2 = @(x) FF2(x) + an*cos(i*pi*x/p);
end

figure(4)
fplot(@(x) f(x),[-p,p],'-b','linewidth',2)
hold on
fplot(@(x) FF1(x),[-p,p],'-.g','linewidth',2)
fplot(@(x) FF2(x),[-p,p],'-or','linewidth',2)
hold off
grid on
title('Example 3: Cosine Expansion','fontsize',14,'fontweight','bold');
xlabel('X','fontsize',12,'fontweight','bold');
ylabel('F(x)','fontsize',12,'fontweight','bold');
set(gca,'fontsize',12,'fontweight','bold');
legend('f(x)','N=3','N=5','location','best');

%% Example 3 - Sine Series Expansion
clear;

f = @(x) x.^2;
L = 2; %<--to plot we need to define this
p = L;

FF1 = @(x) 0;
N1 = 3;

for i = 1:N1
    bn = (2/p)*integral(@(x) f(x).*sin(i*pi*x/p),0,p); %notice bounds of integration.
    
    FF1 = @(x) FF1(x) + bn*sin(i*pi*x/p);
end
FF2 = @(x) FF1(x);
N2 = 5;
for i = (N1+1):N2    
    bn = (2/p)*integral(@(x) f(x).*sin(i*pi*x/p),0,p);
    
    FF2 = @(x) FF2(x) + bn*sin(i*pi*x/p);
end

figure(5)
fplot(@(x)ex3_odd(x),[-p,p],'-b','linewidth',2)
hold on
fplot(@(x) FF1(x),[-p,p],'-.g','linewidth',2)
fplot(@(x) FF2(x),[-p,p],'-or','linewidth',2)
hold off
grid on
title('Example 3: Sine Series','fontsize',14,'fontweight','bold');
xlabel('X','fontsize',12,'fontweight','bold');
ylabel('F(x)','fontsize',12,'fontweight','bold');
set(gca,'fontsize',12,'fontweight','bold');
legend('f(x)','N=5','N=15','location','best');

%% Example 3 - Full Fourier Series
clear
f = @(x) ex3_full(x);

N1 = 5;
p = 2;

a0 = (1/p)*integral(@(x) f(x),-p,p);
FF1 = @(x) a0/2;

for i = 1:N1
    an = (1/p)*integral(@(x) f(x).*cos(i*pi*x/p),-p,p);
    bn = (1/p)*integral(@(x) f(x).*sin(i*pi*x/p),-p,p);
    
    FF1 = @(x) FF1(x) + an*cos(i*pi*x/p) + bn*sin(i*pi*x/p);
end
FF2 = @(x) FF1(x);
N2 = 15;
for i = (N1+1):N2
    an = (1/p)*integral(@(x) f(x).*cos(i*pi*x/p),-p,p);
    bn = (1/p)*integral(@(x) f(x).*sin(i*pi*x/p),-p,p);
    
    FF2 = @(x) FF2(x) + an*cos(i*pi*x/p) + bn*sin(i*pi*x/p);
end

figure(6)
fplot(@(x)f(x),[-p,p],'-b','linewidth',2)
hold on
fplot(@(x) FF1(x),[-p,p],'-.g','linewidth',2)
fplot(@(x) FF2(x),[-p,p],'-or','linewidth',2)
hold off
grid on
title('Example 3: Full Fourier','fontsize',14,'fontweight','bold');
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

function y = ex3_odd(x)
[m,n] = size(x); % write your function to expect vector inputs
y = nan(m,n); % construct y such that it has the same shape as x
for i = 1:length(x)
    if (x(i) > -2) && (x(i) < 0)
        y(i) = -(x(i).^2);
    elseif(x(i)>=0) && (x(i)<2)
        y(i) = x(i).^2;
    end
end

end

function y = ex3_full(x)
[m,n] = size(x); % write your function to expect vector inputs
y = nan(m,n); % construct y such that it has the same shape as x
for i = 1:length(x)
    if (x(i) > -2) && (x(i) < 0)
        y(i) = ((x(i)+2).^2);
    elseif(x(i)>=0) && (x(i)<2)
        y(i) = x(i).^2;
    end
end


end