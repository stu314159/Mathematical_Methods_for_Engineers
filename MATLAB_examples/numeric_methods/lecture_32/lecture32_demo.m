%% Lecture 32 Demo
clear
clc
close 'all'

%% Define constants
a = 0; b = 0.1;
N = 2000;
Ac = 1.6e-5; % m^2, fin cross sectional area
P = 0.016; % m, perimeter of pin cross section
h_c = 40; % W/m^2-K, convective heat transfer coefficient of air around pin
k = 250; % W/m-K, thermal conductivity of pin material
emiss = 0.5; % emissivity of pin material
sigma_sb = 5.67e-8; % W/m^2-K^4, Stefan-Boltzmann constant
Ts = 293; % K, temperature of surrounding air

Ta = 473;
Tb = 293;

%% Discretize the space and get Operators
x = linspace(a,b,N); x = x';

Dx_op = Dx(a,b,N);
Dxx_op = Dxx(a,b,N);

alpha_1 = h_c*P/(k*Ac);
alpha_2 = emiss*sigma_sb*P/(k*Ac);

L = Dxx_op - alpha_1*speye(N,N);

% apply boundary condition:
L(1,:) = 0; L(1,1) = 1; 
L(N,:) = 0; L(N,N) = 1;

%% Estimate initial temperature distribution
To = linspace(Ta,Tb,N); To = To';

%% Form RHS

b = -ones(N,1)*Ts*alpha_1 - ones(N,1)*(Ts^4)*alpha_2;
phi = (To.^4)*alpha_2;

RHS = b - phi;
% apply BC to RHS
RHS(1) = Ta; RHS(N) = Tb;

% Ready to apply iteration
tol = 1e-7;
imax = 100;

T = To;

for i = 1:imax
    
    %solve the system of equations for new estimate of temperature
    Tnew = L\RHS; 
    
    % obtain convergence criterion "error estimate"
    Err = norm(T - Tnew,Inf)/norm(T,Inf);
    
    % exit loop if error is within tolerance
    if Err <tol
        fprintf('Success!! Iteration converged after %d iterations.\n',i);
        fprintf('Error estimate: %g \n',Err);
        break;
    end
    
    % error tolerance not met, prepare for next iteration
    phi = (Tnew.^4)*alpha_2;
    RHS = b - phi;
    RHS(1) = Ta; RHS(N) = Tb; %re-apply BC   
    T = Tnew;
    
end

if i == imax
    fprintf('Error! Solution not converged after %i iterations.\n',imax);
    fprintf('Last residual = %g \n',Err);
end

fprintf('\n\n Plotting Last Solution: \n\n');
figure(1)
plot(x,T,'-r','linewidth',3);
title('Solution of Non-linear BVP with Finite Difference Method',...
    'fontsize',16,'fontweight','bold');
xlabel('X (m) ','fontsize',14,'fontweight','bold');
ylabel('T (K) ','fontsize',14,'fontweight','bold');
set(gca,'fontsize',12,'fontweight','bold');
grid on

%% Local Functions
function dx_sp = Dx(a,b,N)
% function dx_sp = Dx(a,b,N) returns a sparse matrix
% for a first order differentiation matrix using 2nd-order
% centered-difference for interior nodes and 2nd-order
% forward/backward-difference nodes for the respective endpoints of the
% domain.
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
% forward/backward-difference nodes for the respective endpoints of the
% domain.
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
