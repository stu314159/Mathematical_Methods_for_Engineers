%% Double Fourier Series Demonstration
clear
clc
close 'all'

%% Parameters
a = 4;
b = 3;

N = 20;

ex_select = 2;
% 1 smooth
% 2 not smooth

switch ex_select
    case 1
        f = @(x,y) x.*(a-x).*y.*(b-y);
        
    case 2
       f = @(x,y) ex1(x,y,a,b);
       
    otherwise
        error('Invalid Example Choice!');
end

Nx = 250;
Ny = 250;

X = linspace(0,a,Nx);
Y = linspace(0,b,Ny);

[XX,YY] = meshgrid(X,Y);

figure(1)
surf(XX,YY,f(XX,YY),'edgecolor','none');
title('f(x,y)','fontsize',16,'fontweight','bold');
xlabel('X','fontsize',14,'fontweight','bold');
ylabel('Y','fontsize',14,'fontweight','bold');
zlabel('f(X,Y)','fontsize',14,'fontweight','bold');
set(gca,'fontsize',12,'fontweight','bold');

%% Double Fourier Expansion
Xn = @(x,n) sin(n.*pi.*x./a);
Yn = @(y,n) sin(n.*pi.*y./b);

A = nan(N,N);
FF = @(x,y) 0;

for m = 1:N
    for n = 1:N
        A(m,n) = integral2(@(x,y) f(x,y).*Xn(x,m).*Yn(y,n),...
            0,a,0,b)./...
            integral2(@(x,y) (Xn(x,m).^2).*(Yn(y,n).^2),...
            0,a,0,b); %(same as formula - 4/(a*b) above)
        FF = @(x,y) FF(x,y)+A(m,n).*Xn(x,m).*Yn(y,n);
    end
end

figure(2)
surf(XX,YY,FF(XX,YY),'edgecolor','none');
title('Double Fourier Expansion of f(x,y)','fontsize',16,'fontweight','bold');
xlabel('X','fontsize',14,'fontweight','bold');
ylabel('Y','fontsize',14,'fontweight','bold');
zlabel('FF(X,Y)','fontsize',14,'fontweight','bold');
set(gca,'fontsize',12,'fontweight','bold');

%% Local functions
function z = ex1(x,y,a,b)
[mx,nx] = size(x); % write your function to expect vector inputs
[my,ny] = size(y);
% for this implementation, I will expect x and y to have the same size
assert((mx==my) && (nx == ny),...
    'error: x and y must have same size');
z = nan(mx,nx); % construct z to be the same as X and Y

% there are fancier ways to do this, but we will use a very simple
% implementation
for i = 1:mx
    for j = 1:nx
        if (x(i,j) < a/2) && (y(i,j) < b/2)
            z(i,j) = -1;
        else
            z(i,j) = 1;
        end
    end
end
end


