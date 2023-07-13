%% Double Fourier Series Demonstration
clear
clc
close 'all'

%% Parameters
a = 4;
b = 3;
k = 0.3;

N = 10;

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



%% Heat Equation Solution
Xn = @(x,n) sin(n.*pi.*x./a);
Yn = @(y,n) sin(n.*pi.*y./b);
Tn = @(t,m,n) exp(-k.*((m*pi/a).^2 + (n*pi/b).^2).*t);

A = nan(N,N); % note that there is no requirement that n == m
U = @(x,y,t) 0;

for m = 1:N
    for n = 1:N
        % compute the m,n Fourier Coefficient
        A(m,n) = integral2(@(x,y) f(x,y).*Xn(x,m).*Yn(y,n),...
            0,a,0,b)./...
            integral2(@(x,y) (Xn(x,m).^2).*(Yn(y,n).^2),...
            0,a,0,b); 
        % add the term to the expansion
        U = @(x,y,t) U(x,y,t) + A(m,n).*Xn(x,m).*Yn(y,n).*Tn(t,m,n);
    end
end


Nx = 250;
Ny = 250;

X = linspace(0,a,Nx);
Y = linspace(0,b,Ny);

[XX,YY] = meshgrid(X,Y);

Tmax = 1;
NT = 10;
T = linspace(0,Tmax,NT);

figure(1)
for i = 1:NT
    surf(XX,YY,U(XX,YY,T(i)),'edgecolor','none');
    title_str = sprintf('U(x,y,t), t = %g',T(i));
    title(title_str,'fontsize',16,'fontweight','bold');
    xlabel('X','fontsize',14,'fontweight','bold');
    ylabel('Y','fontsize',14,'fontweight','bold');
    zlabel('U(X,Y,T)','fontsize',14,'fontweight','bold');
    set(gca,'fontsize',12,'fontweight','bold');
    
    switch ex_select
        case 1
            axis([0 a 0 b -1 10]);
        case 2
            axis([0 a 0 b -2 2]);            
    end
    
    pause(Tmax/(NT-1));
    
end

%% Local functions
function z = ex1(x,y,a,b)
[mx,nx] = size(x); % write your function to expect vector inputs
[my,ny] = size(y);
% for this implementation, I will expect x and y to have the same size
assert((mx==my) && (nx == ny),'error: x and y must have same size');
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



