%% Preconditioning Example
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

imax = 2000;
tol = 1e-7;
% without preconditioning
[x_jac1,norm_res1,num_iter1,exit_code1] = ...
    jacobi_solver(A,b,x_in,tol,imax);
if exit_code1 == 1
   fprintf('Un-preconditioned Jacobi solution successful! \n');
   fprintf('Number of iterations: %d \n',num_iter1);
   fprintf('tol = %g \n',norm_res1);
end

% suppose we use the LU factorization of A and inverted it.
[L,U] = lu(A);
nnzL = nnz(L); nnzU = nnz(U);

PA = U\(L\A); Pb = U\(L\b);% wouldn't *actually* do this but...

[x_jac2,norm_res2,num_iter2,exit_code2] = ...
    jacobi_solver(PA,Pb,x_in,tol,imax);
if exit_code2 == 1
   fprintf('Preconditioned Jacobi solution successful! \n');
   fprintf('Number of iterations: %d \n',num_iter2);
   fprintf('tol = %g \n',norm_res2);
end


opts_ilu.type='ilutp';
opts_ilu.droptol=1e-1;
[iL,iU] = ilu(A,opts_ilu);
nnz_iL = nnz(iL); nnz_iU = nnz(iU);

iPA = iU\(iL\A); iPb = iU\(iL\b);%wouldn't *actually* do this but...

[x_jac3,norm_res3,num_iter3,exit_code3] = ...
    jacobi_solver(iPA,iPb,x_in,tol,imax);
if exit_code3 == 1
   fprintf('Preconditioned Jacobi solution successful! \n');
   fprintf('Number of iterations: %d \n',num_iter3);
   fprintf('tol = %g \n',norm_res3);
end

%% Local Functions
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

