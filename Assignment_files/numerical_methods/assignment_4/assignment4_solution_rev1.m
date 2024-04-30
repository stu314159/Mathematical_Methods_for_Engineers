%% Assignment #4 Solution
clear
clc
close 'all'

%% Problem #3

fprintf('\n\n Problem #3 \n\n');
A = [4 0 1 0 1;
    2 5 -1 1 0;
    1 0 3 -1 0;
    0 1 0 4 -2;
    1 0 -1 0 5];
b = [32 19 14 -2 41]';

[L,U,P] = lu(A);
bp = P*b;
y = ForwardSub(L,bp);
x_p = BackwardSub(U,y);
fprintf('a) solution: \n'); disp(x_p);

rr = norm(A*x_p - b,2)/norm(b,2);
fprintf('b) relative residual: %g \n',rr);

Ainv = inv(A);

fprintf('c) 2-norm of A: %g \n',norm(A,2));
fprintf('2-norm of inv(A): %g \n',norm(Ainv,2));

eb_lower = 1/(norm(A,2)*norm(Ainv,2))*rr;
eb_upper = norm(Ainv,2)*norm(A,2)*rr;
fprintf('d)\n');
fprintf('lower error bound: %g \n',eb_lower);
fprintf('upper error bound: %g \n',eb_upper);

fprintf('Repeating for nos5.mtx\n');
[A,~,cols,~] = mmread('nos5.mtx');
A = full(A);
% fprintf('Test matrix with %d rows, %d cols \n',...
%     rows,cols);
b = ones(cols,1);
[L,U,P] = lu(A);
bp = P*b;
y = ForwardSub(L,bp);
x_p = BackwardSub(U,y);

rr = norm(A*x_p - b,2)/norm(b,2);
Ainv = inv(A);
eb_lower = 1/(norm(A,2)*norm(Ainv,2))*rr;
eb_upper = norm(Ainv,2)*norm(A,2)*rr;

fprintf('lower error bound: %g \n',eb_lower);
fprintf('upper error bound: %g \n',eb_upper);

%% Problem #4
fprintf('\n\n Problem #4 \n\n');
A = [7, 3, -1, 2;
    3, 8, 1, -4;
    -1, 1, 4, -1;
    2, -4, -1, 6];
b = [-1; 0; -3; 1];

imax = 20000;
tol = 1e-9;
x_in = zeros(4,1);

[x_jac,~,ni_jac,~] = ...
    jacobi_solver(A,b,x_in,tol,imax);
fprintf('Jacobi Solver iterations: %g \n',ni_jac);

[x_gs,~,ni_gs,~] = gs_solver(A,b,x_in,tol,imax);
fprintf('Gauss-Seidel iterations: %g \n',ni_gs);

omega = 1.4;
[x_sor,~,ni_sor,~] = ...
    sor_solver(A,b,x_in,tol,imax,omega);
fprintf('SOR iterations: %g \n',ni_sor);

%% Problem #5
fprintf('\n\n Problem #5 \n\n');
[A,rows,cols,entries] = mmread('nos5.mtx');
%A = full(A);
% fprintf('Test matrix with %d rows, %d cols \n',...
%     rows,cols);
b = ones(cols,1);
x = zeros(cols,1);
tol = 1e-4;
imax = 20000;
%%
min_w = 1.5;
max_w = 1.96;
N_sor = 30;
w = linspace(min_w,max_w,N_sor);
iter_sor = nan(N_sor,1);
for k = 1:N_sor
    [~,~,iter_sor(k),~] = sor_solver(A,b,x,tol,imax,w(k));
    %fprintf('Estimated Relative Error (SOR) = %g \n',ere_sor1);
    fprintf('Number of iterations (SOR) = %g \n',iter_sor(k));
end

figure(1)
plot(w,iter_sor,'-k','linewidth',3)
title('Iterations Required vs. Relaxation Parameter',...
    'FontSize',16,'FontWeight','bold');
grid on
xlabel('Relaxation Parameter [\omega]','fontsize',14,'FontWeight','bold');
ylabel('Iterations','FontSize',14,'FontWeight','bold');
set(gca,'fontsize',12,'fontweight','bold');


% Solve using Built-in Iterative Methods

%% without preconditioning
[x0,fl0,rr0,it0,rv0] = pcg(A,b,tol,imax);

if (fl0 == 0) 
    fprintf('pcg without preconditioner solution successful!\n');
    fprintf('Residual norm: %g, after %d iterations.\n',...
        rr0,it0);
end

if (fl0 == 1)
    fprintf('pcg without preconditionging failed to converge. \n');
end


%% Generalized Minimum residual

[x3,fl3,rr3,it3,rv3] = gmres(A,b,[],tol,min(imax,size(A,1)));
if fl3 == 0
    fprintf('GMRES without preconditioner solution successful!\n');
    fprintf('Residual norm: %g, after %d outer iterations\n',...
        rr3,it3(1));
    fprintf('and %d inner iterations.\n',it3(2));
end

%% add incomplete lu factorization preconditioner
opts_ilu.type='ilutp';
opts_ilu.droptol=1e-4;
% opts_ilu.type='crout';
% opts_ilu.milu = 'row';
% opts_ilu.droptol = 0.01;

[L,U] = ilu(A,opts_ilu); 
[x4,fl4,rr4,it4,rv4] = ...
    gmres(A,b,[],tol,min(imax,size(A,1)),L,U);

if fl4 == 0
    fprintf('GMRES with ilu preconditioner solution successful!\n');
    fprintf('Residual norm: %g, after %d outer iterations\n',...
        rr4,it4(1));
    fprintf('and %d inner iterations.\n',it4(2));
end

%% built-in solve:

f_bi = @() A\b;

time_direct = timeit(f_bi);
fprintf('Using backslash: %g sec \n',time_direct);

%% gmres

f_gmres = @() gmres(A,b,[],tol,min(imax,size(A,1)),L,U);

time_gmres = timeit(f_gmres);
fprintf('Using gmres: %g sec\n',time_gmres);

%% Local Functions
function y = ForwardSub(L,b)
% based on Example 4-6, program 4-4 in the text pg 126-127
n = length(b);
y(1,1) = b(1)/L(1,1);
for i = 2:n
   y(i,1)=(b(i)-L(i,1:(i-1))*y(1:(i-1),1))./L(i,i); 
end
end

function x = BackwardSub(U,y)
% based on Example 4-6, program 4-5 in the text pg 127
n = length(y);
x(n,1) = y(n)/U(n,n);
for i = (n-1):-1:1
   x(i,1)=(y(i)-U(i,(i+1):n)*x((i+1):n,1))/U(i,i); 
end
end


function [x,norm_res,num_iter,exit_code] = ...
    jacobi_solver(A,b,x_in,tol,imax)
[rows,~] = size(A);
rel_update = inf;
x = x_in;
for iter = 1:imax
    for i = 1:rows
       x(i) = (b(i) - A(i,:)*x_in)/A(i,i);
       x(i) = x(i) + x_in(i); %add back j=i term
    end  
    
    if norm(x_in,2) ~= 0 % prevent nan
        % estimated relative error.
        rel_update = norm(x - x_in,2)/norm(x_in,2);
    end
    
    % check exit criteria; if applicable set answer and exit function
    if rel_update < tol
        exit_code = 1; % success
        break; % "break out" of the for loop
    end
    
    if iter == imax
        exit_code = 0; % maximum iterations reached
    end
    x_in = x;
    
end
norm_res = rel_update;
num_iter = iter;
end

function [x,norm_res,num_iter,exit_code] = ...
    gs_solver(A,b,x_in,tol,imax)
[rows,~] = size(A);
rel_update = inf;
x = x_in;
for iter = 1:imax
    for i = 1:rows
       x(i) = b(i)/A(i,i);
       x(i) = x(i) - A(i,1:(i-1))*x(1:(i-1))/A(i,i); % use updated terms
       x(i) = x(i) - A(i,(i+1):end)*x_in((i+1):end)/A(i,i); % use old terms      
    end  
    
    if norm(x_in,2) ~= 0 % prevent nan
        rel_update = norm(x - x_in,2)/norm(x_in,2);
    end
    
    % check exit criteria; if applicable set answer and exit function
    if rel_update < tol
        exit_code = 1; % success
        break; % "break out" of the for loop
    end
    
    if iter == imax
        exit_code = 0; % maximum iterations reached
    end
    x_in = x;
    
end
norm_res = rel_update;
num_iter = iter;
end

function [x,norm_res,num_iter,exit_code] = ...
    sor_solver(A,b,x_in,tol,imax,omega)
[rows,~] = size(A);
rel_update = inf;
x = x_in;
for iter = 1:imax
    for i = 1:rows
       
       x(i) = (1-omega)*x_in(i); 
       x(i) = x(i) + omega*b(i)/A(i,i);
       x(i) = x(i) - omega*A(i,1:(i-1))*x(1:(i-1))/A(i,i); % use updated terms
       x(i) = x(i) - omega*A(i,(i+1):end)*x_in((i+1):end)/A(i,i); % use old terms      
    end  
    
    if norm(x_in,2) ~= 0 % prevent nan
        rel_update = norm(x - x_in,2)/norm(x_in,2);
    end
    
    if (iter > 5) && isnan(rel_update)
       fprintf('SOR nan detected, exiting iteration\n');
       num_iter = nan;
       return;
    end
    
    % check exit criteria; if applicable set answer and exit function
    if rel_update < tol
        exit_code = 1; % success
        break; % "break out" of the for loop
    end
    
    if iter == imax
        exit_code = 0; % maximum iterations reached
    end
    x_in = x;
    
end
norm_res = rel_update;
num_iter = iter;
end
