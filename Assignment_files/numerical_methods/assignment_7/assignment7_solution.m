%% Assignment7_solution.m
clear
clc
close 'all'

%% Problem #1 - 8.17
fprintf('\n\n Problem #1 \n\n');
x1 = [1.1 1.2 1.3 1.4 1.5];
y1 = [0.6133 0.7822 0.9716 1.1814 1.4117];
yd = FirstDeriv(x1,y1);

fprintf('The first derivative at x = %g is %g \n',...
    x1(3),yd(3));

%% Problem #2 - 8.18
fprintf('\n\n Problem #2 \n\n');
x2 = [-1 -0.5 0 0.5 1 1.5 2 2.5 3 3.5 4 4.5];
y2 = [-3.632 -0.3935 1 0.6487 -1.282 -4.518 -8.611 ...
    -12.82 -15.91 -15.88 -9.402 9.017];
ydd = SecDeriv(x2,y2);
fprintf('The second derivative at x = %g is %g \n',...
    x2(8),ydd(8));


%% Problem 3 - 8.37
fprintf('\n\n Problem #3 \n\n');
t = [0 4 8 12 16 20 24 28 ...
    32 36 40 44 48 52 56 60];
r = 1e3*[18.803 18.861 18.946 19.042 ...
    19.148 19.260 19.376 19.495 ...
    19.617 19.741 19.865 19.990 ...
    20.115 20.239 20.362 20.484];
theta = [0.7854 0.7792 0.7701 0.7594 ...
    0.7477 0.7350 0.7215 0.7073 0.6925 ...
    0.6771 0.6612 0.6448 0.6280 0.6107 ...
    0.5931 0.5750];

drdt = FirstDeriv(t,r); 
ddrdt = SecDeriv(t,r);
dtheta_dt = FirstDeriv(t,theta);
ddtheta_dt = SecDeriv(t,theta);

v = sqrt(drdt.^2 + (r.*dtheta_dt).^2);
a1 = ddrdt-r.*(dtheta_dt).^2;
a2 = r.*ddtheta_dt + 2*drdt.*dtheta_dt;
a = sqrt(a1.^2 + a2.^2);

figure(2)
yyaxis left
plot(t,v,'linewidth',3);
ylabel('Velocity [m/s]','fontsize',14,'fontweight','bold');
title('Velocity and Acceleration','fontsize',16,...
    'fontweight','bold');
xlabel('Time [s]','fontsize',14,'fontweight','bold');
set(gca,'fontsize',10,'fontweight','bold');

yyaxis right
plot(t,a,'linewidth',3);
ylabel('Acceleration [m/s^2]','fontsize',14,...
    'fontweight','bold');
grid on
set(gca,'fontsize',10,'fontweight','bold');

%% Problem #4 - 9.19
fprintf('\n\n Problem #4 \n\n');
z = [-18 -12 -6 0 6 12 18];% in, elevation
d = [28 30.2 31.5 32 31.5 30.2 28]; % in, diameter
% use diameter data from solution not from the book's problem statement

S = 2*pi*IntPointsTrap(z,d./2);
V = pi*IntPointsTrap(z,(d./2).^2);

fprintf('Surface Area: %g in^2\n',S);
fprintf('Volume: %g in^3\n',V);

%% Problem #5 - 9.20
fprintf('\n\n Problem #5 \n\n');
x = [0 0.3 0.6 0.9 1.2 1.5 1.8];
f_x = [0.5 0.6 0.8 1.3 2 3.2 4.8];

I = SimpsonPoints(x,f_x);

fprintf('The integral is: %g \n',I);

%% Problem #6 - 9.25
fprintf('\n\n Problem #6 \n\n');
Fun = @(x) exp(-x.^2);
a = 0; b = 3;

GQ5ab = @(Fun1,a1,b1) GaussQuad1D(Fun1,a1,b1,5);

I6 = GaussQuad5ab(Fun,a,b);
fprintf('The integral is: %g \n',I6);

I6a = GaussQuad1D(Fun,a,b,5);
fprintf('Using GaussQuad1D: %g \n',I6a);

I6b = GQ5ab(Fun,a,b);
fprintf('Using GQ5ab: %g \n',I6b);


%% Local Function
function yd = FirstDeriv(x,y)
assert(min(size(x))==1);% x and y should both be vectors
assert(min(size(y))==1);
assert(max(size(x))==max(size(y))); % x and y the same length
assert(max(size(y))>=3); % enough data for proposed method
[m,n] = size(y);
yd = nan(m,n); % construct output same shape as y.
h = x(2)-x(1); % it would be good to verify this
for s = 2:(length(h)-1)
    ht = x(s+1)-x(s);
    assert(abs(h - ht)<=eps);
end
% but might be going overboard.

% left end
yd(1) = (-3*y(1) + 4*y(2) - y(3))/(2*h);

% right end
yd(end) = (y(end-2) - 4*y(end-1) + 3*y(end))/(2*h);

% middle
i = 2:(length(x)-1);
yd(i) = (y(i+1) - y(i-1))/(2*h);

end

function ydd = SecDeriv(x,y)
assert(min(size(x))==1);% x and y should both be vectors
assert(min(size(y))==1);
assert(max(size(x))==max(size(y))); % x and y the same length
assert(max(size(y))>=3); % enough data for proposed method
[m,n] = size(y);
ydd = nan(m,n); % construct output same shape as y.
h = x(2)-x(1); % it would be good to verify this
for s = 2:(length(h)-1)
    ht = x(s+1)-x(s);
    assert(abs(h - ht)<=eps);
end
% but might be going overboard.

% 4-point forward difference on the left end point
ydd(1) = (2*y(1)-5*y(2)+4*y(3)-y(4))/(h*h);

% 4-point backward difference on the right end point
ydd(end) = (-y(end-3)+4*y(end-2)-5*y(end-1)+2*y(end))/(h*h);

% 3-point centered difference everywhere in between
i = 2:(length(y)-1);
ydd(i) = (y(i-1)-2*y(i)+y(i+1))/(h*h);

end

function I = IntPointsTrap(x,y)
assert(min(size(x))==1);% x and y should both be vectors
assert(min(size(y))==1);
assert(max(size(x))==max(size(y))); % x and y the same length
h = x(2:end) - x(1:(end-1));
fip = y(2:end); fi = y(1:(end-1));
I = sum((h/2).*(fip + fi));
end

function I = SimpsonPoints(x,y)
assert(min(size(x))==1);% x and y should both be vectors
assert(min(size(y))==1);
assert(max(size(x))==max(size(y))); % x and y the same length
assert(mod(length(x),2)==1); % vector has odd number of points
assert(length(x)>=3); % at least 3 points given.
h = x(2) - x(1);
if (length(x) > 3)
    % general composite Simpson's Rule
    I = (h/3)*(y(1)+y(end)+sum(4*y(2:2:end))...
        +sum(2*y(3:2:(end-1))));
else
    % special case for only 3 points
    I = (h/3)*(y(1)+4*y(2)+y(3));
end
end

function I = GaussQuad5ab(Fun,a,b)
% check a > b and fix it if you need to.
if (a > b)
    tmp = b;
    b = a;
    a = tmp;
end

% set sample points and weights for 5-point GQ.
t = nan(5,1);
t(1) = -0.90617985;
t(5) = -t(1);
t(2) = -0.53846931;
t(4) = -t(2);
t(3) = 0;

w = nan(5,1);
w(1) = 0.2369269;
w(5) = w(1);
w(2) = 0.4786287;
w(4) = w(2);
w(3) = 0.5688889;

% set up mapping from [-1,1] -> [a,b]
xT = @(t) ((b-a)*t + a + b)/2;
Jac = (b - a)/2;

x = xT(t);

I = w'*Fun(x)*Jac;
end

function [intF, xgl, wgl] = GaussQuad1D(F,a,b,P)
%myGaussQuad1D(f,a,b,P) performs P-point Gaussian Quadrature of function F
%over interval [a,b]
%   input:  F - function handle
%           a - lower limit of integration
%           b - upper limit of integration
%           P - # of Gauss Points to use in integration
%  output: intF - numeric integral of F
%           xgl - vector of gauss points
%           wgl - vector of weights
%
%
%% Generate Legendre Polynomials of Order 0 through P
% Store handles to these functions in a cell array.
isQuadSet = false; % flag for special exit of quadrature scheme.
Pn = cell(P+1,1);
Pn{1} = @(x) 1;
Pn{2} = @(x) x;
if P == 0
    error('P must be greater than 0');
elseif P == 1
    xgl = 0; % for P == 1, GQ reduces to midpoint rule over entire domain.
    wgl = 2;
    isQuadSet = true;
else
    for n = 2:P
        % use recurrence relation to generate higher order Legendre
        % Polynomials ("Pn functions")
        Pn{n+1} = @(x) (2*(n-1)+1)*x.*Pn{n}(x)./((n-1)+1) - (n-1)*Pn{n-1}(x)./((n-1)+1);
    end
end

if ~(isQuadSet)
    %% Compute Roots to the Pth order Legendre Polynomial
    
    % get an approximate value of the zeros from the Chebychev points
    Tch = @(n) cos(((2*(1:n)) - 1)*pi./(2*n));
    xEst = Tch(P);
    
    % use fzero and approximate root to find root of the Pn polynomial.
    xgl = NaN(1,P);
    for r = 1:P
        xgl(r) = fzero(Pn{P+1},xEst(r));
    end
    
    %% Sample lower order Pn functions at roots of Pth order function
    % These values form the matrix that will be used to find the weights for
    % the quadrature method.
    % form coefficient matrix
    if P == 1
        A = xgl(1);
    else
        A = NaN(P,P);
        A(1,:) = Pn{1}(xgl);
        A(2,:) = Pn{2}(xgl);
        for n = 2:(P-1)
            A((n+1),:) = Pn{n+1}(xgl);
        end
    end
    
    %% Form LHS vector
    % These are equal to the integral of the lower order Pn functions over the
    % domain.  For P0, the integral equals 2; for all other orders, the
    % integral is zero.
    k = zeros(P,1); k(1) = 2;
    
    %% Solve for the weights
    wgl = A\k;
end

%% Change Variables to Scale Interval to [-1,1]
xT = @(t) ((b-a)*t + a + b)/2;
Jac = (b - a)/2;

%% Perform the Integration
intF = F(xT(xgl))*wgl*Jac;

end




