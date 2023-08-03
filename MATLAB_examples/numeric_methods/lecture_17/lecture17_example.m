%% Lecture 17 examples
clear
clc
close 'all'

%% Example 1
f = @(x) cos(x);
df = @(x) -sin(x);
n = 7;
Nx = 2^n;
xMin = 0; xMax = pi;
X = linspace(xMin,xMax,Nx);
h = X(2) - X(1);

df_numeric = nan(1,Nx);
% use fwd difference at left end
df_numeric(1) = (1/h)*(f(X(2))- f(X(1)));
% use backward difference at right end
df_numeric(end) = (1/h)*(f(X(end)) - f(X(end-1)));
% use centered difference everywhere else
i = 2:(Nx-1); ip = i+1; im = i-1;
df_numeric(2:(end-1)) = (1/(2*h))*(f(X(ip))-f(X(im)));

% plot the function and its derivative
Xp = linspace(xMin,xMax,10000);
figure(1)
plot(Xp,f(Xp),'-b',...
    X,df_numeric,'-r','linewidth',2);
grid on;
title('First Derivative of $\cos{x}$',...
    'fontsize',14,...
    'fontweight','bold','Interpreter','latex');
xlabel('X','fontsize',12,'fontweight','bold');
ylabel('Y','fontsize',12,'fontweight','bold');
legend('f(x)','f^{\prime}(x)');
set(gca,'fontsize',12,'fontweight','bold');

%% Get Convergence Rate
N = 5:15;

rel_err = nan(1,length(N));
h_err = nan(1,length(N));

for s = 1:length(N)
    Nx = 2^N(s);
    xMin = 0; xMax = pi;
    X = linspace(xMin,xMax,Nx);
    h = X(2) - X(1);
    h_err(s) = h;
    
    df_numeric = nan(1,Nx);
    % use fwd difference at left end
    df_numeric(1) = (1/h)*(f(X(2))- f(X(1)));
    % use backward difference at right end
    df_numeric(end) = (1/h)*(f(X(end)) - f(X(end-1)));
    % use centered difference everywhere else
    i = 2:(Nx-1); ip = i+1; im = i-1;
    df_numeric(i) = (1/(2*h))*(f(X(ip))-f(X(im)));
    
    % estimate the relative error
    x_err = 1:Nx; % include the end points
    %x_err = 2:(Nx-1); % exclude the end points
   
    rel_err(s) = norm(df_numeric(x_err) - df(X(x_err)),2)...
        /norm(df(X(x_err)),2);    
end

% add gauge lines for convergence
c1 = 0;
h1 = h_err + c1;

c2 = 0;
h2 = h_err.^2 + c2;


figure(2)
loglog(h_err,rel_err,'-b',...
    h_err,h1,'--r',...
    h_err,h2,'--g','linewidth',3);
title('Convergence Behavior','fontsize',14,...
    'fontweight','bold');
xlabel('h','fontsize',12,'fontweight','bold');
ylabel('Relative Error','fontsize',12,'fontweight','bold');
grid on
set(gca,'fontsize',10,'fontweight','bold');
legend('Estimate','h^1 convergence','h^2 convergence',...
    'location','best');

%% Repeat example
% this time using 3-point fwd and bwd difference operators

N = 5:15;

rel_err = nan(1,length(N));
h_err = nan(1,length(N));

for s = 1:length(N)
    Nx = 2^N(s);
    xMin = 0; xMax = pi;
    X = linspace(xMin,xMax,Nx);
    h = X(2) - X(1);
    h_err(s) = h;
    
    df_numeric = nan(1,Nx);
    % use fwd difference at left end
    df_numeric(1) = (1/(2*h))*(-3*f(X(1))+4*f(X(2))-f(X(3)));
    % use backward difference at right end
    df_numeric(end) = (1/(2*h))*(3*f(X(end))-4*f(X(end-1))+...
        f(X(end-2)));
    % use centered difference everywhere else
    i = 2:(Nx-1); ip = i+1; im = i-1;
    df_numeric(i) = (1/(2*h))*(f(X(ip))-f(X(im)));
    
    % estimate the relative error
    x_err = 1:Nx; % include the end points
    %x_err = 2:(Nx-1); % exclude the end points
   
    rel_err(s) = norm(df_numeric(x_err) - df(X(x_err)),2)...
        /norm(df(X(x_err)),2);    
end

% add gauge lines for convergence
c1 = 0;
h1 = h_err + c1;

c2 = 0;
h2 = h_err.^2 + c2;


figure(3)
loglog(h_err,rel_err,'-b',...
    h_err,h1,'--r',...
    h_err,h2,'--g','linewidth',3);
title('Convergence Behavior','fontsize',14,...
    'fontweight','bold');
xlabel('h','fontsize',12,'fontweight','bold');
ylabel('Relative Error','fontsize',12,'fontweight','bold');
grid on
set(gca,'fontsize',10,'fontweight','bold');
legend('Estimate','h^1 convergence','h^2 convergence',...
    'location','best');

%% Matrix Representation 

N = 10;
Nx = 2^N;
xMin = 0; xMax = pi;
X = linspace(xMin,xMax,Nx);
h = X(2) - X(1);

A = zeros(Nx,Nx); % matrix for 2nd derivative operator

% set coefficients for first equation
A(1,1) = 2/h^2; A(1,2) = -5/h^2;
A(1,3) = 4/h^2; A(1,4) = -1/h^2;

% set coefficents for all interior equations
for m=2:(Nx-1)
    A(m,m) = -2/h^2; A(m,m-1) = 1/h^2;
    A(m,m+1) = 1/h^2;
end

% set coefficents for last equation
A(Nx,Nx) = 2/h^2; A(Nx,Nx-1) = -5/h^2;
A(Nx,Nx-2) = 4/h^2; A(Nx,Nx-3) = -1/h^2;

% get finite difference second derivative 
fd_ddf = A*(f(X)'); % note syntax 

ddf = @(x) -cos(x);


% make a plot
figure(4)
plot(Xp,ddf(Xp),'-r',...
    X,fd_ddf,'--b',...
    'linewidth',2);
title('Numeric 2nd Derivative','fontsize',14,...
    'fontweight','bold');
xlabel('X','fontsize',12,'fontweight','bold');
ylabel('Y','fontsize',12,'fontweight','bold');
grid on
legend('Analytic','Numeric','location','best')


%% Example 3 Boundary Value Problem Example

N = 10;
Nx = 2^N;
xMin = 0; xMax = pi;
X = linspace(xMin,xMax,Nx);
h = X(2) - X(1);

A = zeros(Nx,Nx); % matrix for 2nd derivative operator

% set coefficients for first equation
A(1,1) = 2/h^2; A(1,2) = -5/h^2;
A(1,3) = 4/h^2; A(1,4) = -1/h^2;

% set coefficents for all interior equations
for m=2:(Nx-1)
    A(m,m) = -2/h^2; A(m,m-1) = 1/h^2;
    A(m,m+1) = 1/h^2;
end

% set coefficents for last equation
A(Nx,Nx) = 2/h^2; A(Nx,Nx-1) = -5/h^2;
A(Nx,Nx-2) = 4/h^2; A(Nx,Nx-3) = -1/h^2;

% Set boundary conditions on first and last equation
A(1,:) = 0; A(1,1) = 1;
A(end,:) = 0; A(end,end) = 1;


%  use MATLAB tools to get eigenvalues and eigenvectors of A
[V,D] = eigs(A);


figure(5)
subplot(4,1,1)
plot(X,(V(:,1))/norm(V(:,1),2),'.b')

subplot(4,1,2)
plot(X,(V(:,2))/norm(V(:,2),2),'.b')

subplot(4,1,3)
plot(X,(V(:,3))/norm(V(:,3),2),'.b')

subplot(4,1,4)
plot(X,(V(:,4))/norm(V(:,4),2),'.b')
