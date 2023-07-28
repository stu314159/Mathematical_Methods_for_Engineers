%% Example: Jacobi Method Demonstration
clear
clc
close 'all'

sys_choice = 2;

switch sys_choice
    
    case 1
        A = [9 -2 3 2;
            2 8 -2 3;
            -3 2 11 -4;
            -2 3 2 10];
        b = [54.5;
            -14;
            12.5;
            -21];
        rows = 4; cols = 4;
        x_in = zeros(cols,1);
        show_x = 1;
           
    case 2
        [A,rows,cols,entries] = mmread('gr_30_30.mtx');
        fprintf('Test matrix with %d rows, %d cols \n',...
            rows,cols);
        b = rand(cols,1);
        x_in = zeros(rows,1);
        figure
        spy(A)
        title('Sparsity Pattern of Test Matrix')
        show_x = 0;

    otherwise
        error('Invalid system choice!\n');
end

%% Solve using Jacobi Iteration

% set stopping criteria
imax = 1500; % max iterations
tol = 1e-7; % relative size of update

fprintf('Calculation with simple Jacobi iteration \n');
tic;
[x_jac,norm_res,num_iter,exit_code] = ...
    jacobi_solver(A,b,x_in,tol,imax);
time_jac = toc;

if exit_code == 1
    fprintf('Iterative solution successful!\n');
    fprintf('Residual norm: %g, after %d iterations.\n',...
        norm_res,num_iter);
end

if show_x
   fprintf('x = \n'); disp(x_jac); 
end

fprintf('Calculation with vectorized Jacobi iteration \n');
tic;
[x_jac_v,norm_res_v,num_iter_v,exit_code_v] = ...
    jacobi_solver_v(A,b,x_in,tol,imax);
time_jacv = toc;
if exit_code_v == 1
    fprintf('Iterative solution converged!\n');
    fprintf('Residual norm: %g, after %d iterations.\n',...
        norm_res_v,num_iter_v);
end

if show_x
   fprintf('x = \n'); disp(x_jac_v); 
end


%% Check Relative Residual

fprintf('Relative residual Jacobi: %g \n',norm_res);
fprintf('Relative residual vectorized Jacobi: %g \n',...
    norm_res_v);


%% Compare Relative Error
x_gold = A\b;

rel_error_jacobi = norm(x_jac_v-x_gold,2)/norm(x_gold,2);
fprintf('Relative error: %g \n',...
    rel_error_jacobi);

%% Compare execution time

fprintf('\n');
fprintf('Time for simple Jacobi: %g sec\n',time_jac);
fprintf('Time for vectorized Jacobi: %g sec\n',time_jacv);

%% Local functions
function [x_new,norm_res,num_iter,exit_code] = ...
    jacobi_solver(A,b,x_in,tol,imax)
[n,~] = size(A);
rel_update = inf;

x_new = x_in; % initialize x_new
for iter = 1:imax
    if (iter > 1) && (mod(iter,500) == 0)
        fprintf('Iteration: %d, relative update:%g. \n',...
            iter,rel_update);
    end
 
    for i = 1:n
        x_new(i) = (1/A(i,i))*(b(i) - ...
            A(i,1:(i-1))*x_in(1:(i-1)) - ...
            A(i,(i+1):n)*x_in((i+1):n));
    end

    if norm(x_in,"inf") ~= 0 % prevent nan
        rel_update = ...
            norm(x_new - x_in,"inf")/norm(x_in,"inf");
    end
    
    % check exit criteria
    if rel_update < tol
        exit_code = 1; % success
        break; % "break out" of the for loop
    end
    
    if iter == imax
        % maximum iterations reached
        exit_code = 0; 
    end
    x_in = x_new;
    
end
norm_res = norm(A*x_new - b,2)/norm(x_new,2);
num_iter = iter;
end


function [x_new,norm_res,num_iter,exit_code] = ...
    jacobi_solver_v(A,b,x_in,tol,imax)
rel_update = inf;
x_new = x_in; % initialize x_new
K = -(A - diag(diag(A)));
M_inv = sparse(diag(1./diag(A)));

for iter = 1:imax
    if (iter > 1) && (mod(iter,500) == 0)
        fprintf('Iteration: %d, relative update:%g. \n',...
            iter,rel_update);
    end

    x_new = M_inv*(K*x_in + b);

    if norm(x_in,"inf") ~= 0 % prevent nan
        rel_update = ...
            norm(x_new - x_in,"inf")/norm(x_in,"inf");
    end
    
    % check exit criteria
    if rel_update < tol
        exit_code = 1; % success
        break; % "break out" of the for loop
    end
    
    if iter == imax
        % maximum iterations reached
        exit_code = 0; 
    end
    x_in = x_new;
    
end
norm_res = norm(A*x_new - b,2)/norm(x_new,2);
num_iter = iter;
end
