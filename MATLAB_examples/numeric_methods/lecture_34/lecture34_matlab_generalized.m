%% Finite Element Example in 1D
clear
clc
close 'all'

%% Exact Solution
u_exact = @(x) x - sinh(x)./sinh(1);
a = 0; b = 1;
x_gold = linspace(a,b,1000);

Ua = 0; Ub = 0; % Boundary Conditions
%% Finite Element Parameters and Data Structures
nelem = 3; % select # of elements
order = 1;

[gcoord,nodes] = genMesh1D(a,b,nelem,order);

% nldofs = number of local dofs for each element
nldofs = order+1;

% global x-coordinate of each node
%gcoord = linspace(a,b,nelem+1); % x-coordinates of all nodes
nnodes = length(gcoord);

% sample points for Gauss Quadrature
nqp = order+1; % number of sample points requested
[q,w] = getGPandWeights(nqp);

% initialize global arrays
K1 = zeros(nnodes,nnodes);
K2 = zeros(nnodes,nnodes);
R = zeros(nnodes,1);

% carry out assembly process
for ele = 1:nelem

    % local arrays to be populated
    k1 = zeros(nldofs,nldofs);
    k2 = zeros(nldofs,nldofs);
    r = zeros(nldofs,1);
    
    % local mapping for GQ
    aL = gcoord(nodes(ele,1));
    bL = gcoord(nodes(ele,end));
    xT = @(t) ((bL - aL)*t + aL + bL)/2;
    Jac = (bL - aL)/2;

    % Get sample points for shape functions
    xgl = gcoord(nodes(ele,:));
    
    % get Lagrange Interpolant of requested order
    H = getLagrangeInterp(xgl);
    Hp = getLagrangeInterpDeriv(xgl);

    for qp = 1:nqp
        
        % sum weighted contribution at Gauss Points
        for i = 1:nldofs
            for j = 1:nldofs
                k1(i,j) = k1(i,j) + Hp{i}(xT(q(qp)))*Hp{j}(xT(q(qp)))*w(qp);
                k2(i,j) = k2(i,j) + H{i}(xT(q(qp)))*H{j}(xT(q(qp)))*w(qp);
            end % j
            r(i) = r(i) + xT(q(qp))*H{i}(xT(q(qp)))*w(qp);
        end % i
     
    end % qp
    % apply Jacobian to map to physical coordinates
    k1 = k1*Jac;
    k2 = k2*Jac;
    r = r*Jac;

    % add local arrays to global arrays "assembly"
    for i = 1:nldofs
        for j = 1:nldofs
            row = nodes(ele,i); col = nodes(ele,j);
            K1(row,col) = K1(row,col) + k1(i,j);
            K2(row,col) = K2(row,col) + k2(i,j);
        end % j
        dof = nodes(ele,i);
        R(dof) = R(dof) + r(i);
    end % i

end % numel

% apply boundary conditions
A = -K1 - K2;
RHS = -R;
A(1,:) = 0; A(1,1) = 1; RHS(1) = Ua;
A(nnodes,:) = 0; A(nnodes,nnodes) = 1; RHS(nnodes) = Ub;

% solve the system of equations
u = A\RHS;

%% Plot the solution and check error
x = gcoord;
figure(1)
subplot(2,1,1)
plot(x,u','-b',...
    x_gold,u_exact(x_gold),'--r',...
    'linewidth',2);
title('Finite Element Solution',...
    'FontSize',14,'FontWeight','bold');
%xlabel('X','FontSize',12,'FontWeight','bold');
ylabel('U(X)','FontSize',12,'FontWeight','bold');
grid on
legend('U - FEM','U - Exact','Location','best');

err = u' - u_exact(x);
subplot(2,1,2)
plot(x,err,'-b','linewidth',2);
xlabel('X','FontSize',12,'FontWeight','bold');
grid on

rel_err = norm(err,2)/norm(u_exact(x),2);
fprintf('Relative error: %g \n',rel_err);


%% Local Functions
function [xgl, wgl] = getGPandWeights(P)
%GetGP(P) Get Gauss Points and Weights for P-point Gauss Quad
%   input:  P - # of Gauss Points to use in integration
%  output:  xgl - vector of gauss points
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
end

function H = getLagrangeInterp(Xi)
%function dH = getLagrangeInterp(Xi)
% input Xi - vector of sample points
% output H - cell array containing Lagrange Interp Functions

n = length(Xi);
H = cell(n,1);

for i = 1:n
    L = @(x) 1; %<-- initialize the Lagrange Function
    for j = 1:n
        if j ~= i
            L = @(x) L(x).*((x - Xi(j))./(Xi(i) - Xi(j)));
        end
    end
    H{i} = L;
end

end


function dH = getLagrangeInterpDeriv(Xi)
%function dH = getLagrangeInterpDerive(Xi)
% input: Xi - vector of sample points
% output dH - cell array containing the derivative of Lagrange Interpolant
% functions.
n = length(Xi);
dH = cell(n,1);

for i = 1:n
    dLi = @(x) 0;
    for k = 1:n
        if k ~= i
            dLp = @(x) 1;
            for j = 1:n
                if ((j ~= i) && (j ~= k))
                    dLp = @(x) dLp(x).* (x - Xi(j))./(Xi(i) - Xi(j));
                end
            end
            dLi = @(x) dLi(x) + (1./(Xi(i) - Xi(k))).*dLp(x);
        end
        
    end
    dH{i} = dLi;
end

end

function [L0,L0_1,L0_2] = legendre_poly(p,x)
%---------------------------------------------------------------------%
%This code computes the Legendre Polynomials and its 1st and 2nd derivatives
%Written by F.X. Giraldo on 4/2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
L1=0;L1_1=0;L1_2=0;
L0=1;L0_1=0;L0_2=0;

for i=1:p
   L2=L1;L2_1=L1_1;L2_2=L1_2;
   L1=L0;L1_1=L0_1;L1_2=L0_2;
   a=(2*i-1)/i;
   b=(i-1)/i;
   L0=a*x*L1 - b*L2;
   L0_1=a*(L1 + x*L1_1) - b*L2_1;
   L0_2=a*(2*L1_1 + x*L1_2) - b*L2_2;
end
end


function [xgl,wgl] = legendre_gauss_lobatto(P)
%---------------------------------------------------------------------%
%This code computes the Legendre-Gauss-Lobatto points and weights
%which are the roots of the Lobatto Polynomials.
%Written by F.X. Giraldo on 4/2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%

p=P-1; %Order of the Polynomials
ph=floor( (p+1)/2 );
xgl = nan(1,P);
wgl = nan(1,P);
for i=1:ph
   % estimate the roots
   x=cos( (2*i-1)*pi/(2*p+1) );
   for k=1:20
      % evaluate L, dL, ddL
      [L0,L0_1,L0_2]=legendre_poly(p,x); %Compute Nth order Derivatives of Legendre Polys
      
      %Get new Newton Iteration
      dx=-(1-x^2)*L0_1/(-2*x*L0_1 + (1-x^2)*L0_2);
      x=x+dx;
      if (abs(dx) < 1.0e-20) 
         break
      end
   end
   xgl(p+2-i)=x;
   wgl(p+2-i)=2/(p*(p+1)*L0^2);
end

%Check for Zero Root
if (p+1 ~= 2*ph)
   x=0;
   [L0,~,~]=legendre_poly(p,x);
   xgl(ph+1)=x;
   wgl(ph+1)=2/(p*(p+1)*L0^2);
end
   
%Find remainder of roots via symmetry
for i=1:ph
   xgl(i)=-xgl(p+2-i);
   wgl(i)=+wgl(p+2-i);
end
   
end

function [gcoord,nodes] = genMesh1D(a,b,nelem,order)
%genMesh1D(a,b,nelem,order) generates the global coordinates of all mesh
%points and provides a mapping between global node number and local
%element dof for all elements.

ldofs = order + 1; % for 1st order need 2 points, 2nd order 3 pts, etc...

ndofs = nelem+1 + (ldofs-2)*nelem;
gcoord = nan(1,ndofs);
nodes = nan(nelem,ldofs);
ele_boundaries = linspace(a,b,nelem+1);
dof_it = 1;
for e = 1:nelem 
    aL = ele_boundaries(e);
    bL = ele_boundaries(e+1);
    xT = @(t) ((bL - aL)*t + aL + bL)/2;
    [xgl,~]=legendre_gauss_lobatto(ldofs);
    xgl = xT(xgl); % translate to current element

    for ld = 1:ldofs
        nodes(e,ld) = dof_it;
        gcoord(dof_it) = xgl(ld);
        dof_it = dof_it + 1;
    end % ld
    dof_it = dof_it - 1; % set back for shared node at element boundary
end % e

end

