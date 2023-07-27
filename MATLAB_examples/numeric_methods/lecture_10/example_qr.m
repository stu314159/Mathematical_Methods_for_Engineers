%% Example QR Decomposition
clear
clc
close 'all'

sys_choice = 4;

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
        error('Invalid system choice!\n');
end

[Q,R] = Example_QR(A);
x = BackwardSub(R,Q'*b);

if show_x
   fprintf('Q: \n');
   disp(Q);
   fprintf('R: \n');
   disp(R);
   fprintf('QR: \n');
   disp(Q*R);
   fprintf('A: \n');
   disp(A);
   
   fprintf('x: \n');
   disp(x);
   
end

%% Compare Relative Error
x_gold = A\b;

rel_error_qr = norm(x-x_gold,2)/norm(x_gold,2);
fprintf('Relative error QR: %g \n',...
    rel_error_qr);


%% Local functions
function [Q,R] = Example_QR(A)
% QR factorization based on Modified Graham-Schmidt
% Implementation based on Trefethen Algorithm 8.1

[m,n] = size(A); % get rows and columns of A
R = zeros(m,n);
Q = zeros(m,m);
V = A;
for k = 1:n
    R(k,k) = norm(V(:,k),2); % magnitude of column k
    Q(:,k) = V(:,k)./R(k,k); % kth column of Q is normalized
    for j = (k+1):n
        % get fraction of column k in direction of column j
        R(k,j) = Q(:,k)'*V(:,j); 
        % subtract that fraction out from column j
        V(:,j) = V(:,j) - R(k,j)*Q(:,k);        
    end        
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
