%% Finite Element Example in 1D
clear
clc
close 'all'

%% Exact Solution
u_exact = @(x) x - sinh(x)./sinh(1);
a = 0; b = 1;
x_gold = linspace(a,b,1000);

Ua = 0; Ub = 0;
%% Finite Element Parameters and Data Structures
nelem = 3; % select # of elements

% global x-coordinate of each node
gcoord = linspace(a,b,nelem+1); % x-coordinates of all nodes
nnodes = length(gcoord);

% array to hold node numbers for each
% node in every element
nodes = nan(nelem,2);
% assign node numbers to each element
nodes(:,1) = 1:nelem;% x_i node number
nodes(:,2) = 2:(nelem+1);%x_i+1 node number
% this acts like a map between local node numbers for each element
% and global node numbers that tie to the overall geometry

% sample points for Gauss Quadrature
q = [-0.57735027; 0.57735027];
% weights
w = [1; 1];
nqp = length(q); % number of quadrature points

% initialize global arrays
K1 = zeros(nnodes,nnodes);
K2 = zeros(nnodes,nnodes);
R = zeros(nnodes,1);

% carry out assembly process
for ele = 1:nelem

    % local arrays to be populated
    k1 = zeros(2,2);
    k2 = zeros(2,2);
    r = zeros(2,1);
    
    % local mapping for GQ
    aL = gcoord(nodes(ele,1));
    bL = gcoord(nodes(ele,2));
    xT = @(t) ((bL - aL)*t + aL + bL)/2;
    Jac = (bL - aL)/2;

    % shape functions for this element
    H = cell(2,1);
    hi = bL - aL;
    H{1} = @(x) (bL-x)/hi;
    H{2} = @(x) (x-aL)/hi;

    Hp = cell(2,1);
    % for generality, use functions
    Hp{1} = @(x) -1/hi;
    Hp{2} = @(x) 1/hi;


    for qp = 1:nqp
        
        % sum weighted contribution at Gauss Points
        k1(1,1) = k1(1,1) + Hp{1}(xT(q(qp)))*Hp{1}(xT(q(qp)))*w(qp);
        k1(1,2) = k1(1,2) + Hp{1}(xT(q(qp)))*Hp{2}(xT(q(qp)))*w(qp);
        k1(2,1) = k1(2,1) + Hp{2}(xT(q(qp)))*Hp{1}(xT(q(qp)))*w(qp);
        k1(2,2) = k1(2,2) + Hp{2}(xT(q(qp)))*Hp{2}(xT(q(qp)))*w(qp);
        
        k2(1,1) = k2(1,1) + H{1}(xT(q(qp)))*H{1}(xT(q(qp)))*w(qp);
        k2(1,2) = k2(1,2) + H{1}(xT(q(qp)))*H{2}(xT(q(qp)))*w(qp);
        k2(2,1) = k2(2,1) + H{2}(xT(q(qp)))*H{1}(xT(q(qp)))*w(qp);
        k2(2,2) = k2(2,2) + H{2}(xT(q(qp)))*H{2}(xT(q(qp)))*w(qp);
        
        r(1) = r(1) + xT(q(qp))*H{1}(xT(q(qp)))*w(qp);
        r(2) = r(2) + xT(q(qp))*H{2}(xT(q(qp)))*w(qp);        

    end % qp
    % apply Jacobian to map to physical coordinates
    k1 = k1*Jac;
    k2 = k2*Jac;
    r = r*Jac;

    % add local arrays to global arrays ("assembly")
    for i = 1:2
        for j = 1:2
            row = nodes(ele,i); col = nodes(ele,j);
            K1(row,col) = K1(row,col) + k1(i,j);
            K2(row,col) = K2(row,col) + k2(i,j);
        end % j
        dof = nodes(ele,i);
        R(dof) = R(dof) + r(i);
    end % i

end % numel

% Gather into system matrix and vector
A = -K1 - K2;
RHS = -R;

% apply boundary conditions
A(1,:) = 0; A(1,1) = 1; RHS(1) = Ua;
A(nnodes,:) = 0; A(nnodes,nnodes) = 1; RHS(nnodes) = Ub;

% solve the system of equations
u = A\RHS;

%% Plot the solution and check error
x = gcoord;
figure(1)
plot(x,u','-b',...
    x_gold,u_exact(x_gold),'--r',...
    'linewidth',2);
title('Finite Element Solution',...
    'FontSize',14,'FontWeight','bold');
xlabel('X','FontSize',12,'FontWeight','bold');
ylabel('U(X)','FontSize',12,'FontWeight','bold');
grid on
legend('U - FEM','U - Exact','Location','best');




