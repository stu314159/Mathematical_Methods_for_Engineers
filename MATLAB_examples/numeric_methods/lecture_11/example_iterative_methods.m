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
        error('Invalide system choice!\n');
end




%% Solve using simple Jacobi Iteration
imax = 4;
tol = 1e-4;

tic;
[x_jac,norm_res,num_iter,exit_code] = ...
    jacobi_solver_simple(A,b,x_in,tol,imax);
t1 = toc;

tic;
[~,~,~,~] = ...
    jacobi_solver(A,b,x_in,tol,imax);
t2 = toc;

if show_x
   fprintf('x = \n'); disp(x_jac); 
end

fprintf('Time with non-vectorized solver: %g \n',t1);
fprintf('time with vectorized solver: %g \n',t2);

%% Now just use Vectorized Solver
% compare Jacobi, Gauss Seidel, and SOR
imax = 2000;
tol = 1e-7;


[x_jac2,norm_res2,num_iter2,exit_code2] = ...
    jacobi_solver(A,b,x_in,tol,imax);

%% Use Vectorized Gauss Seidel Solver

[x_gs,nr_gs,it_gs,ec_gs] = ...
    gs_solver(A,b,x_in,tol,imax);

%% Use SOR Solver
omega = 1.85;
[x_sor,nr_sor,it_sor,ec_sor] = ...
    sor_solver(A,b,x_in,tol,imax,omega);

fprintf('Number of iterations for \n');
fprintf('Jacobi: %d \n',num_iter2);
fprintf('Gauss Seidel: %d \n',it_gs);
fprintf('SOR: %d \n',it_sor);
%% Compare Relative Error
x_gold = A\b;

rel_error_jacobi = norm(x_jac2-x_gold,2)/norm(x_gold,2);
fprintf('Relative error: %g \n',...
    rel_error_jacobi);


%% Local functions
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
    jacobi_solver_simple(A,b,x_in,tol,imax)
[rows,~] = size(A);
rel_update = inf;
x = x_in;
for iter = 1:imax
    for i = 1:rows
        x(i) = b(i)/A(i,i);
        for j = 1:(i-1)
           x(i) = x(i) - (A(i,j)/A(i,i))*x_in(j); 
           
        end
        for j = (i+1):rows
           x(i) = x(i) - (A(i,j)/A(i,i))*x_in(j); 
           
        end
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

