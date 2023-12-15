%% Lecture 20 demo
clear
clc
close 'all'

f = @(x) exp(-x.^2);
a = -3; b = 3;

N = 16;

% you can  pick and choose which outputs you want
% to keep.
intF = GaussQuad1D(f,a,b,N);
[intF1,xgl] = GaussQuad1D(f,a,b,N);
[intF2,xgl2,wgl2] = GaussQuad1D(f,a,b,N);
[~,~,wgl3] = GaussQuad1D(f,a,b,N);


%% Local Functions
function [intF, xgl, wgl] = GaussQuad1D(F,a,b,P)
%GaussQuad1D(f,a,b,P) performs P-point Gaussian quadrature 
%of function F over interval [a,b]
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
    xgl = 0; % for P == 1, GQ reduces to midpoint rule 
    wgl = 2;
    isQuadSet = true;
else
    for n = 2:P
        % use recurrence relation to generate higher 
        % order Legendre Polynomials ("Pn functions")
        Pn{n+1} = @(x) ...
            (2*(n-1)+1)*x.*Pn{n}(x)./((n-1)+1) ...
            - (n-1)*Pn{n-1}(x)./((n-1)+1);
    end
end

if ~(isQuadSet)
    %% Compute Roots to the Pth order 
    % Legendre Polynomial
    
    % get an approximate value of the zeros 
    % from the Chebychev points
    Tch = @(n) cos(((2*(1:n)) - 1)*pi./(2*n));
    xEst = Tch(P);
    
    % use fzero and approximate root to find root of 
    % the Pn polynomial.
    xgl = NaN(1,P);
    for r = 1:P
        xgl(r) = fzero(Pn{P+1},xEst(r));
    end
    
    %% Sample lower order Pn functions at 
    % roots of Pth order function
    % These values form the matrix that will be 
    % used to find the weights for
    % the quadrature method.
    % form coefficient matrix
    if P == 1
        A = xgl(1);
    else
        A = NaN(P,P);
        % A(1,:) = Pn{1}(xgl);
        % A(2,:) = Pn{2}(xgl);
        % for n = 2:(P-1)
        %     A((n+1),:) = Pn{n+1}(xgl);
        % end
        for n = 0:(P-1)
            A((n+1),:) = Pn{n+1}(xgl);
        end
    end
    
    %% Form LHS vector
    % These are equal to the integral of the lower 
    % order Pn functions over the domain.  
    % For P0, the integral equals 2; for all other 
    % orders, the integral is zero.
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

