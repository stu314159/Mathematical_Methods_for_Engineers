%% Lecture 26 MATLAB - Laplace's Eqn
clear
clc
close 'all'

%% Set Parameters
a = 5;
b = 3;
N = 15;

fx_pick = 1;
%[1 | 2]
switch fx_pick
    case 1
        f = @(x) ex1(x,a);
    case 2
        f = @(x) ex2(x,a);
    otherwise
        error('Invalid case!');        
end

%% define the eigenvalues and eigenfunctions
alpha = @(n) n.*pi./a;
F = @(x,n) cos(alpha(n).*x);
G = @(y,n) sinh(alpha(n).*y);

%% Compute coefficients
% compute Ao
Ao = (1/(a*b))*integral(@(x) f(x),0,a);

% initialize solution
u = @(x,y) Ao.*y;
for n = 1:N
    % compute An
    An = (2./(a*G(b,n))).*...
        integral(@(x) f(x).*F(x,n),0,a);
    % update the approximate solution
    u = @(x,y) u(x,y) + An*F(x,n).*G(y,n);
end

%% Make discrete spatial coordinate axes
Nx = ceil(100*a);% "ceil" rounds up to next highest integer
Ny = ceil(100*b);
X = linspace(0,a,Nx);
Y = linspace(0,b,Ny);

[XX,YY] = meshgrid(X,Y);

%% Plot the solution in a 2D plot using surf
figure(1)
surf(XX,YY,u(XX,YY),'edgecolor','none');
title("Lecture 26 Laplace's Equation Example",...
    'fontsize',16,'fontweight','bold');
xlabel('X','fontsize',14,'fontweight','bold');
ylabel('Y','fontsize',14,'fontweight','bold');
zlabel('u(X,Y)','fontsize',14,'fontweight','bold');
set(gca,'fontsize',12,'fontweight','bold');

%% Local functions
function y = ex1(x,a)
[m,n] = size(x);
y = nan(m,n);
for i = 1:length(x)
    if(x(i)>= 0) && (x(i)<a/2)
        y(i) = x(i).^2;
    elseif(x(i) >= a/2) && (x(i)<a)
        y(i) = (a/2).^2;
    end
end    
end

function y = ex2(x,a)
[m,n] = size(x);
y = nan(m,n);
for i = 1:length(x)
    if(x(i)>= 0) && (x(i)<a/2)
        y(i) = x(i);
    else
        y(i) = a - x(i);
    end
end    
end