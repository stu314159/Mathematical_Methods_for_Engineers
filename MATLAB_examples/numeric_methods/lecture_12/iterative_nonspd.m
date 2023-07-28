%% Example: Iterative Methods
clear
clc
close 'all'

sys_choice = 1;

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
        [A,rows,cols,entries] = mmread('e05r0000.mtx');
        % est condition number 6e4
        fprintf('Test matrix with %d rows, %d cols \n',...
            rows,cols);
        b = rand(cols,1);
        x_in = zeros(rows,1);
        figure
        spy(A)
        title('Sparsity Pattern of Test Matrix')
        show_x = 0;
        
    case 3
        [A,rows,cols,entries] = mmread('e20r0000.mtx');
        % est condition number 5.45e10
        fprintf('Test matrix with %d rows, %d cols \n',...
            rows,cols);
        b = rand(cols,1);
        x_in = zeros(rows,1);
        figure
        spy(A)
        title('Sparsity Pattern of Test Matrix')
        show_x = 0;
        
    case 4
        [A,rows,cols,entries] = mmread('impcol_d.mtx');
        % est condition number 1.9e3
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
%% Solve using Built-in Iterative Methods
imax = 20000;
tol = 1e-4;



%% Generalized Minimum residual

[x3,fl3,rr3,it3,rv3] = gmres(A,b,[],tol,min(imax,size(A,1)));
if fl3 == 0
    fprintf('GMRES without preconditioner solution successful!\n');
    fprintf('Residual norm: %g, after %d outer iterations\n',...
        rr3,it3(1));
    fprintf('and %d inner iterations.\n',it3(2));
end

% add incomplete lu factorization preconditioner
opts_ilu.type='ilutp';
opts_ilu.droptol=1e-4;

[L,U] = ilu(A,opts_ilu); 
[x4,fl4,rr4,it4,rv4] = ...
    gmres(A,b,[],tol,min(imax,size(A,1)),L,U);

if fl4 == 0
    fprintf('GMRES with ilu preconditioner solution successful!\n');
    fprintf('Residual norm: %g, after %d outer iterations\n',...
        rr4,it4(1));
    fprintf('and %d inner iterations.\n',it4(2));
end



%% Local functions
% function [x,norm_res,num_iter,exit_code] = ...
%     jacobi_solver(A,b,x_in,tol,imax)
% [rows,~] = size(A);
% rel_update = inf;
% x = x_in;
% for iter = 1:imax
%     if (iter > 1) && (mod(iter,1000) == 0)
%         fprintf('Iteration: %d, relative update:%g. \n',...
%             iter,rel_update);
%     end
%     for i = 1:rows
%        x(i) = (b(i) - A(i,:)*x_in)/A(i,i);
%        x(i) = x(i) + x_in(i); %add back j=i term
%     end  
%     
%     if norm(x_in,2) ~= 0 % prevent nan
%         rel_update = norm(x - x_in,2)/norm(x_in,2);
%     end
%     
%     % check exit criteria; if applicable set answer and exit function
%     if rel_update < tol
%         exit_code = 1; % success
%         break; % "break out" of the for loop
%     end
%     
%     if iter == imax
%         exit_code = 0; % maximum iterations reached
%     end
%     x_in = x;
%     
% end
% %norm_res = norm(A*x - b,2)/norm(x,2);
% norm_res = rel_update;
% num_iter = iter;
% end

% function [x,norm_res,num_iter,exit_code] = ...
%     sor_solver(A,b,x_in,tol,imax,omega)
% [rows,~] = size(A);
% rel_update = inf;
% x = x_in;
% for iter = 1:imax
%     for i = 1:rows
%        
%        x(i) = (1-omega)*x_in(i); 
%        x(i) = x(i) + omega*b(i)/A(i,i);
%        x(i) = x(i) - omega*A(i,1:(i-1))*x(1:(i-1))/A(i,i); % use updated terms
%        x(i) = x(i) - omega*A(i,(i+1):end)*x_in((i+1):end)/A(i,i); % use old terms      
%     end  
%     
%     if norm(x_in,2) ~= 0 % prevent nan
%         rel_update = norm(x - x_in,2)/norm(x_in,2);
%     end
%     
%     % check exit criteria; if applicable set answer and exit function
%     if rel_update < tol
%         exit_code = 1; % success
%         break; % "break out" of the for loop
%     end
%     
%     if iter == imax
%         exit_code = 0; % maximum iterations reached
%     end
%     x_in = x;
%     
% end
% norm_res = rel_update;
% num_iter = iter;
% end
% 




