%% Example: Iterative Methods
clear
clc
close 'all'

sys_choice = 6;

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
        A = sparse(A);
           
    case 2
        [A,rows,cols,entries] = mmread('nos1.mtx');
        % condition number est 2.5e7
        fprintf('Test matrix with %d rows, %d cols \n',...
            rows,cols);
        b = ones(cols,1);
        figure
        spy(A)
        title('Sparsity Pattern of Test Matrix')
        show_x = 0;
        x_in = zeros(cols,1);
        
    case 3
        [A,rows,cols,entries] = mmread('nos2.mtx');
        % est condition number 6.3e9
        fprintf('Test matrix with %d rows, %d cols \n',...
            rows,cols);
        b = ones(cols,1);
        figure
        spy(A)
        title('Sparsity Pattern of Test Matrix')
        show_x = 0;
        x_in = zeros(cols,1);
        
    case 4
        [A,rows,cols,entries] = mmread('nos7.mtx');
        % condition number 4.1e9
        fprintf('Test matrix with %d rows, %d cols \n',...
            rows,cols);
        b = ones(cols,1);
        figure
        spy(A)
        title('Sparsity Pattern of Test Matrix')
        show_x = 0;
        x_in = zeros(cols,1);
        
    case 5
        [A,rows,cols,entries] = mmread('gr_30_30.mtx');
        % est condition number 3.8e2 with weak diagonal dominance
        fprintf('Test matrix with %d rows, %d cols \n',...
            rows,cols);
        b = ones(cols,1);
        figure
        spy(A)
        title('Sparsity Pattern of Test Matrix')
        show_x = 0;
        x_in = zeros(cols,1);

    case 6
        A = delsq(numgrid('S',102));
        b = ones(size(A,1),1);
        

    otherwise
        error('Invalid system choice!\n');
end




%% Solve using Built-in Iterative Methods
imax = 100;
tol = 1e-8;

% without preconditioning
[x0,fl0,rr0,it0,rv0] = pcg(A,b,tol,imax);

if (fl0 == 0) 
    fprintf('pcg without preconditioner solution successful!\n');
    fprintf('Residual norm: %g, after %d iterations.\n',...
        rr0,it0);
end

if (fl0 == 1)
    fprintf('pcg without preconditionging failed to converge. \n');
end

%% preconditioned conjugate gradient
opts.type = 'ict';
opts.droptol = 1e-4;
opts.michol = 'on';
L = ichol(A,opts);
%L = ichol(A);
[x1,fl1,rr1,it1,rv1] = pcg(A,b,tol,imax,L,L');

if fl1 == 0
    fprintf('pcg with ichol preconditioner solution successful!\n');
    fprintf('Residual norm: %g, after %d iterations.\n',...
        rr1,it1);
end


%% Generalized Minimum residual

[x3,fl3,rr3,it3,rv3] = gmres(A,b,[],tol,min(imax,size(A,1)));
if fl3 == 0
    fprintf('GMRES without preconditioner solution successful!\n');
    fprintf('Residual norm: %g, after %d outer iterations\n',...
        rr3,it3(1));
    fprintf('and %d inner iterations.\n',it3(2));
end

%% GMRES with Incomplete LU Preconditioner
opts_ilu.type='ilutp';
opts_ilu.droptol=1e-3;
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

%% Compare to SOR Iteration
omega = 1.8;
[x2,norm_res,num_iter,exit_code] = ...
    sor_solver(A,b,x_in,tol,imax,omega);

if exit_code == 1
    fprintf('SOR iteration solution successful!\n');
    fprintf('Residual norm: %g, after %d iterations.\n',...
        norm_res,num_iter);
end
%% Local functions
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

