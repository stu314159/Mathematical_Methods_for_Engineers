function intF = myGaussQuad1Dcomposite(F,a,b,P,N)
%myGaussQuad1D(f,a,b,P) performs P-point Gaussian Quadrature of function F
%over interval [a,b]
%   input:  F - function handle 
%           a - lower limit of integration
%           b - upper limit of integration
%           P - # of Gauss Points to use in integration
%           N - # of subintervals to use
% Written by: CDR Stu Blair, USNA Mechanical Engineering Department
% Date: 25 July 2014
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

%% Perform the integration over each subinterval
x_space = linspace(a,b,N+1);

intF = 0;
for k = 1:N
    xT = @(t) ((x_space(k+1)-x_space(k))*t + x_space(k) + x_space(k+1))/2;
    Jac = (x_space(k+1)-x_space(k))/2;
    intF = intF + F(xT(xgl))*wgl*Jac;
end

end

