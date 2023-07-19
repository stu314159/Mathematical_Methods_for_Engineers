%% Example: LU Factorization
clear
clc
close 'all'

sys_choice = 5;

switch sys_choice
    case 1

        A = [4, -2, -3, 6;
            -6, 7, 6.5, -6;
            1, 7.5, 6.25, 5.5;
            -12, 22, 15.5, -1];
        
        b = [12;
            -6.5;
            16;
            17];
        
        show_x = 1;
        
    case 2
        A = rand(50,50);
        b = rand(50,1);
        
        show_x = 0;
        
    case 3
        A = rand(150,150);
        b = rand(150,1);
        
        show_x = 0;
        
    case 4
        [A,rows,cols,entries] = mmread('lns__131.mtx');
        A = full(A); % convert to dense matrix
        fprintf('Test matrix with %d rows, %d cols \n',...
            rows,cols);
        b = rand(cols,1);
        figure
        spy(A)
        title('Sparsity Pattern of Test Matrix')
        show_x = 0;
        
    case 5
        [A,rows,cols,entries] = mmread('lns__511.mtx');
        A = full(A); % convert to dense matrix
        fprintf('Test matrix with %d rows, %d cols \n',...
            rows,cols);
        b = rand(cols,1);
        figure
        spy(A)
        title('Sparsity Pattern of Test Matrix')
        show_x = 0;

    otherwise
        error('Invalide system choice!\n');
end

[L,U] = LU_Factor(A);
y = ForwardSub(L,b);
x = BackwardSub(U,y);

if show_x
   fprintf('L: \n');
   disp(L);
   fprintf('U: \n');
   disp(U);
   fprintf('LU: \n');
   disp(L*U);
   fprintf('A: \n');
   disp(A);
   
   fprintf('x: \n');
   disp(x);
   
end

%% Part 2: LU with Pivoting
[L,U,P] = lu(A); % use MATLAB built-in function
bp = P*b;
y = ForwardSub(L,bp);
x_p = BackwardSub(U,y); 

%% Compare Relative Error
x_gold = A\b;

rel_error_nopivot = norm(x-x_gold,2)/norm(x_gold,2);
fprintf('Relative error with no pivoting: %g \n',...
    rel_error_nopivot);

rel_error_pivot = norm(x_p-x_gold,2)/norm(x_gold,2);
fprintf('Relative error with pivoting: %g \n',...
    rel_error_pivot);



%% Local functions
function [L,U] = LU_Factor(A)
% LU factorization without pivoting
% Implementation based on Trefethen Algorithm 20.1

[m,n] = size(A); % get rows and columns of A
U = A;
L = eye(m,n); % constructor for identity matrix.

for k = 1:(m-1)
    for j = (k+1):m
        L(j,k) = U(j,k)/U(k,k);
        U(j,k:m) = U(j,k:m) - L(j,k)*U(k,k:m);
    end
end
end

function y = ForwardSub(L,b)
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
