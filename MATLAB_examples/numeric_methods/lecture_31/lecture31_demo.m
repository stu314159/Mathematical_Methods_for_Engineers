%% Lecture 31 Demo
clear
clc
close 'all'

%% Test of 1st Order Differential
a = 0;
b = 10;
N = 100;

x = linspace(a,b,N);
x = x'; %<-- get the shape right.

f = @(x) cos(x);
df = @(x) -sin(x);

y = f(x);
dy_exact = df(x);

Dx_op = Dx(a,b,N); %<-- get a differentiation matrix

dy_num = Dx_op*y;

figure(1)
subplot(3,1,1)
plot(x,y,'linewidth',3)
title('Dx Test');
ylabel('f(x)','fontweight','bold');
grid on

subplot(3,1,2)
plot(x,dy_exact,'linewidth',3)
ylabel('f^{\prime}(x) exact','fontweight','bold');
set(gca,'fontweight','bold');
grid on

subplot(3,1,3)
plot(x,dy_num,'linewidth',3)
ylabel('f^{\prime}(x) numerical','fontweight','bold');
set(gca,'fontweight','bold');
grid on

%% Test 2nd-order differential

ddf = @(x) -cos(x);
dyy_exact = ddf(x);

Dxx_op = Dxx(a,b,N);
dyy_num = Dxx_op*y;

figure(2)
subplot(3,1,1)
plot(x,y,'linewidth',3)
title('Dxx Test');
ylabel('f(x)','fontweight','bold');
grid on;

subplot(3,1,2)
plot(x,dyy_exact,'linewidth',3);
ylabel('f^{\prime\prime}(x) exact','fontweight','bold');
set(gca,'fontweight','bold');
grid on

subplot(3,1,3)
plot(x,dyy_num,'linewidth',3);
ylabel('f^{\prime\prime}(x) numerical','fontweight','bold');
set(gca,'fontweight','bold');
grid on;

%% Solve a mildly interesting problem
%%
% 
%  Find the solution to the problem:
% 
% $$y^{\prime \prime} - 4y = 0$$
% 
% 
% 
%  where y(0) = 0 and y(1) = 5;
% 

a = 0; b = 1;
N = 400;

x = linspace(a,b,N); x = x';

% known exact solution
y_exact = @(x) 5*sinh(2*x)./sinh(2);

Dxx_op = Dxx(a,b,N); %<-- get 2nd order differentiation matrix
L = Dxx_op - 4*speye(N,N);%<-- form the differential operator
rhs = zeros(N,1); %<-- initialize the right hand side

% apply boundary conditions
L(1,:) = 0; L(1,1) = 1; rhs(1) = 0; 
L(N,:) = 0; L(N,N) = 1; rhs(N) = 5;

y = L\rhs;

figure(3)
subplot(2,1,1)
xs = x(1:10:end);
plot(xs,y_exact(xs),'sr',...
    x,y,'-c','linewidth',3);
ylabel('Numeric Solution','fontweight','bold');
title('Solution to: y^{\prime\prime} - 4y = 0;  y(0)=0; y(1)=5','fontsize',14);
set(gca,'fontweight','bold');
legend('Exact','Numeric','location','northwest');
grid on

subplot(2,1,2)
plot(x,abs(y - y_exact(x)),'-r','linewidth',3);
ylabel('Numerical Error','fontweight','bold')
set(gca,'fontweight','bold');
grid on

%% Solve a somewhat more difficult problem
%%
% 
% $$x^2 y^{\prime \prime} - 3xy^{\prime} + 3y = 24x^5$$
% 
%%
% 
%  y(1) = y(2) = 0
% 

a = 1; b = 2;
N = 100;

Dxx_op = Dxx(a,b,N);
Dx_op = Dx(a,b,N);

x = linspace(a,b,N); x = x';

% known exact solution
y_exact = @(x) 12*x - 15*(x.^3) + 3*(x.^5);

L = (sparse(diag(x.^2)))*Dxx_op - ...
    3*sparse(diag(x))*Dx_op + ...
    3*speye(N,N);
rhs = 24*(x.^5);

% apply BCs
L(1,:) = 0; L(1,1) = 1; rhs(1) = 0;
L(N,:) = 0; L(N,N) = 1; rhs(N) = 0;

% solve the system
y = L\rhs;

figure(4)
subplot(2,1,1)
xs = x(1:10:end);
plot(xs,y_exact(xs),'sr',...
    x,y,'-c','linewidth',3);
title('Solution to: x^2y^{\prime\prime} - 3xy^{\prime}+3y=24x^5;  y(1)=y(2)=0',...
    'fontweight','bold','fontsize',14);
ylabel('Numeric Solution','fontweight','bold');
set(gca,'fontweight','bold');
legend('Exact','Numeric');
grid on

subplot(2,1,2)
plot(x,abs(y - y_exact(x)),'-r','linewidth',3);
ylabel('Numerical Error')
set(gca,'fontweight','bold');
grid on

%% Local Functions
function dx_sp = Dx(a,b,N)
% function dx_sp = Dx(a,b,N) returns a sparse matrix
% for a first order differentiation matrix using 2nd-order
% centered-difference for interior nodes and 2nd-order
% forward/backward-difference nodes for the respective 
% endpoints of the domain.
%
% Inputs:
% a - scalar, left endpoint of domain
% b - scalar, right endpoint of domain
% N - number of points in the domain inclusive of the endpoints

% compute the number of entries in the sparse matrix:
% 3-each for the 2 endpoints + 2 each for the N-2 interior points
NumEntries = 3*2 + 2*(N-2);

% Initialize the sparse matrix data vectors
dx_row = nan(NumEntries,1);
dx_col = nan(NumEntries,1);
dx_val = nan(NumEntries,1);

h = (b-a)/(N-1);

% first three entries for the left end-point
dx_row(1) = 1; dx_col(1) = 1; dx_val(1) = -3/(2*h);
dx_row(2) = 1; dx_col(2) = 2; dx_val(2) = 4/(2*h);
dx_row(3) = 1; dx_col(3) = 3; dx_val(3) = -1/(2*h);

ind = 4;

for i = 2:(N-1)
   dx_row(ind) = i; dx_col(ind) = i-1; dx_val(ind) = -1/(2*h);
   ind = ind+1;
   dx_row(ind) = i; dx_col(ind) = i+1; dx_val(ind) = 1/(2*h);
   ind = ind+1;   
    
end

% last three entries for the right end-point
dx_row(ind) = N; dx_col(ind) = N; dx_val(ind) = 3/(2*h);
ind = ind+1;
dx_row(ind) = N; dx_col(ind) = N-1; dx_val(ind) = -4/(2*h);
ind = ind+1;
dx_row(ind) = N; dx_col(ind) = N-2; dx_val(ind) = 1/(2*h);

dx_sp = sparse(dx_row,dx_col,dx_val,N,N);

end

function dxx_sp = Dxx(a,b,N)
% function dxx_sp = Dxx(a,b,N) returns a sparse matrix
% for a second order differentiation matrix using 2nd-order
% centered-difference for interior nodes and 2nd-order
% forward/backward-difference nodes for the respective 
% endpoints of the domain.
%
% Inputs:
% a - scalar, left endpoint of domain
% b - scalar, right endpoint of domain
% N - number of points in the domain inclusive of the endpoints

% compute the number of entries in the sparse matrix:
% 4-each for the 2 endpoints + 3 each for the N-2 interior points
NumEntries = 4*2 + 3*(N-2);

% Initialize the sparse matrix data vectors
dx_row = nan(NumEntries,1);
dx_col = nan(NumEntries,1);
dx_val = nan(NumEntries,1);

h = (b-a)/(N-1);

% first three entries for the left end-point
dx_row(1) = 1; dx_col(1) = 1; dx_val(1) = 2/(h^2);
dx_row(2) = 1; dx_col(2) = 2; dx_val(2) = -5/(h^2);
dx_row(3) = 1; dx_col(3) = 3; dx_val(3) = 4/(h^2);
dx_row(4) = 1; dx_col(4) = 4; dx_val(4) = -1/(h^2);

ind = 5;

for i = 2:(N-1)
   dx_row(ind) = i; dx_col(ind) = i-1; dx_val(ind) = 1/(h^2);
   ind = ind+1;
   dx_row(ind) = i; dx_col(ind) = i; dx_val(ind) = -2/(h^2);
   ind = ind+1;
   dx_row(ind) = i; dx_col(ind) = i+1; dx_val(ind) = 1/(h^2);
   ind = ind+1;   
    
end

% last four entries for the right end-point
dx_row(ind) = N; dx_col(ind) = N; dx_val(ind) = 2/(h^2);
ind = ind+1;
dx_row(ind) = N; dx_col(ind) = N-1; dx_val(ind) = -5/(h^2);
ind = ind+1;
dx_row(ind) = N; dx_col(ind) = N-2; dx_val(ind) = 4/(h^2);
ind = ind+1;
dx_row(ind) = N; dx_col(ind) = N-3; dx_val(ind) = -1/(h^2);

dxx_sp = sparse(dx_row,dx_col,dx_val,N,N);

end
