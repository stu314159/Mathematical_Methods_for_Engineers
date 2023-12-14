function intF =  gaussLegendre(F,a,b,P,N)

[xgl,wgl]=gaussLegendreBasis(P);
%% Perform the integration over each subinterval
x_space = linspace(a,b,N+1);

intF = 0;
for k = 1:N
    % map interval to [-1,1]
    xT = @(t) ((x_space(k+1)-x_space(k))*t + x_space(k) + x_space(k+1))/2;
    Jac = (x_space(k+1)-x_space(k))/2;
    % integrate the subinterval.
    intF = intF + F(xT(xgl))'*wgl*Jac;
end

end


%% function [xgl,wgl] = gaussLegendreBasis(P)
% computes sample points and weights for Gauss-Legendre Quadrature.
% This implementation makes maximum use of mathematical identities and
% built-in features of MATLAB.
function [xgl,wgl] = gaussLegendreBasis(P)

%% compute sample points and weights

% sample points are roots of the P-th order Legendre Polynomial
xgl = sort(roots(legendre_coeff(P)));
wgl = NaN(P,1);
for p = 1:P
   % from Wolfram Mathworld Legendre-Gauss Quadrature eq 13
   wgl(p) = 2*(1-xgl(p)^2)/((((P+1)^2)*(polyval(legendre_coeff(P+1),xgl(p)))^2));
end

end

function c = legendre_coeff(P)
% get the coefficients for the P-th order Legendre Polynomial
c = zeros(P+1,1); 
kL = floor(P/2);
for k = 0:kL
    order = P-2*k;
    ind = (P+1)-(order);
    % from Wolfram MathWorld Legendre Polynomial eq 32
    c(ind) = nchoosek(2*P-2*k,P)*nchoosek(P,k)*((-1)^k)*(1/2)^P;
end

end