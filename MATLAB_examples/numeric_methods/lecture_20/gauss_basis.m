% gauss_basis.m
% Written by: CDR Stu Blair, USNA Mechanical Engineering Department
% Date: 25 July 2014
% Purpose: Compute sample points and weights for P-point Gaussian
% Quadrature Rule using Legendre Polynomials as the basis functions.
%

function [xgl,wgl] = gauss_basis(P)
%%function [xgl, wgl] = gauss_basis(P)
%
% input: P - number of Gauss points. Positive integer.
%
% outputs:  xgl - vector of sample points
%           wgl - vector of weight values
%
% Resulting quadrature rule will integrate continuous functions up to Order
% 2P-1 exactly.

%% Generate Legendre Polynomials of Order 0 through P
% Store handles to these functions in a cell array.
Pn = cell(P+1,1);
Pn{1} = @(x) 1;
Pn{2} = @(x) x;
if P == 0
    error('P must be greater than 0');
elseif P == 1
    xgl = 0; % for P == 1, GQ reduces to midpoint rule over entire domain.
    wgl = 2;
    return
else
    for n = 2:P
        % use recurrence relation to generate higher order Legendre
        % Polynomials ("Pn functions")
        Pn{n+1} = @(x) (2*(n-1)+1)*x.*Pn{n}(x)./((n-1)+1) - (n-1)*Pn{n-1}(x)./((n-1)+1);
    end
end

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
% form my coefficient matrix
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
b = zeros(P,1); b(1) = 2;

%% Solve for the weights
wgl = A\b;

end