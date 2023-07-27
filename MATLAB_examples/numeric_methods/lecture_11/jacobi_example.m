%% Example: Iterative Methods
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
        [A,rows,cols,entries] = mmread('bcsstk14.mtx');
        fprintf('Test matrix with %d rows, %d cols \n',...
            rows,cols);
        b = rand(cols,1);
        x_in = zeros(rows,1);
        figure
        spy(A)
        title('Sparsity Pattern of Test Matrix')
        show_x = 0;

    otherwise
        error('Invalide system choice!\n');
end




%% Solve using Jacobi Iteration
imax = 500;
tol = 1e-6;

[x_jac,norm_res,num_iter,exit_code] = ...
    jacobi_solver(A,b,x_in,tol,imax);

if exit_code == 1
    fprintf('Iterative solution successful!\n');
    fprintf('Residual norm: %g, after %d iterations.\n',...
        norm_res,num_iter);
end

if show_x
   fprintf('x = \n'); disp(x_jac); 
end



%% Compare Relative Error
x_gold = A\b;

rel_error_jacobi = norm(x_jac-x_gold,2)/norm(x_gold,2);
fprintf('Relative error: %g \n',...
    rel_error_jacobi);


%% Local functions
function [x,norm_res,num_iter,exit_code] = ...
    jacobi_solver(A,b,x_in,tol,imax)
[rows,~] = size(A);
rel_update = inf;
x = x_in;
for iter = 1:imax
    if (iter > 1) && (mod(iter,20) == 0)
        fprintf('Iteration: %d, relative update:%g. \n',...
            iter,rel_update);
    end
    for i = 1:rows
       x(i) = (b(i) - A(i,:)*x_in)/A(i,i);
       x(i) = x(i) + x_in(i); %add back j=i term
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
norm_res = norm(A*x - b,2)/norm(x,2);
num_iter = iter;
end

