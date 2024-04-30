%% Assignment #3 Solution
clear
clc
close 'all'

%% Problem #1
fprintf('Problem #1 \n\n');

xMin = -10; xMax = 10; Nx = 500;
yMin = -10; yMax = 10; Ny = 500;

f1 = @(x,y) x.^2 + y.^2 - 25;
df1x = @(x,y) 2*x;
df1y = @(x,y) 2*y;

f2 = @(x,y) x.^2 - y - 2;
df2x = @(x,y) 2*x;
df2y = @(x,y) -1;

F = @(x,y) [f1(x,y); f2(x,y)];

Jac = @(x,y) [df1x(x,y) df1y(x,y); df2x(x,y) df2y(x,y)];

Xo = [2.5, 2.5];
tol = 1e-9;
maxIt = 25;

X_out = NewtonSystemSol(F,Jac,Xo,tol,maxIt);

fprintf('Root found at x = %12.11f, y=%12.11f \n',X_out(1),X_out(2));

X = linspace(xMin,xMax,Nx);
Y = linspace(yMin,yMax,Ny);

[XX,YY] = meshgrid(X,Y);

figure(1)
contour(XX,YY,f1(XX,YY),[0 0],'linewidth',3);
hold on
contour(XX,YY,f2(XX,YY),[0 0],'linewidth',3);
hold off
title('Non-Linear System of Equations','fontsize',14,'fontweight','bold')
grid on
xlabel('X','fontsize',12,'fontweight','bold');
ylabel('Y','fontsize',12,'fontweight','bold');
set(gca,'fontsize',10,'fontweight','bold');



%% Problem #2
fprintf('\n\n Problem #2 \n\n');

%X : [Tc, Jc, Th, Jh]
F = @(x) prob3p40(x);

maxIt = 50;
tol = 1e-9;
options = optimoptions('fsolve','Display','iter-detailed',...
    'MaxIterations',maxIt,'StepTolerance',tol);

Xo = [298, 3000, 298, 5000];

[x,fval,exitflag,output] = fsolve(F,Xo,options);

% fprintf('Root found at: \n'); disp(x);
% fprintf('fval = \n'); disp(fval);
fprintf('exitflag = %d \n',exitflag);
fprintf('Returned values: \n');
fprintf('Tc = %11.8f K \n',x(1));
fprintf('Jc = %11.6f   \n',x(2));
fprintf('Th = %11.8f K \n',x(3));
fprintf('Jh = %11.6f   \n',x(4));


%% Problem #4
fprintf('\n\n Problem #4 \n\n');

A = [2 1 -1;
    1 2 1;
    -1 1 -1];
b = [1; 8; -5];

x4 = Gauss(A,b);

fprintf('x = \n'); disp(x4);

%% Problem #5
fprintf('\n\n Problem #5 \n\n');

A = [0.707 1 0 0 0 0 0 0 0 0 0 0 0;
    -0.707 0 1 0 0 0 0 0 0 0 0 0 0;
    0.7071 0 0 1 0 0 0 0 0 0 0 0 0;
    0 -1 0 0 0.659 1 0 0 0 0 0 0 0; % note error in book solution
    0 0 0 -1 -0.753 0 0 0 0 0 0 0 0;
    0 0 -1 0 -0.659 0 1 0 0 0 0 0 0;
    0 0 0 0 0.753 0 0 1 0 0 0 0 0;
    0 0 0 0 0 -1 0 0 0.659 1 0 0 0;
    0 0 0 0 0 0 0 -1 -0.753 0 0 0 0;
    0 0 0 0 0 0 -1 0 -0.659 0 1 0 0;
    0 0 0 0 0 0 0 0 0.753 0 0 1 0;
    0 0 0 0 0 0 0 0 0 -1 0 0 0.707;
    0 0 0 0 0 0 0 0 0 0 0 -1 -0.7071];

b = [0; 2000; -6229; 0; 600; 0; 0; 0; 800; 0; 2429; 0; 600];

x4 = GaussPivot(A,b);
fprintf('x = \n');
format short e
disp(x4);
format short

%% Local Functions
function out = prob3p40(x)
% expectations on x: [Tc, Jc, Th, Jh]

[m,n] = size(x); % expect scalar or vector input
assert(min(m,n)==1,'Error!  Vector input expected for x! \n');
assert(max(m,n)==4,'Error!  Vector of length 4 is required! \n');
out = nan(m,n); % construct output

out(1) = 5.67e-8*x(1).^4 + 17.41*x(1) - x(2) - 5188.18;
out(2) = x(2) - 0.71*x(4) + 7.46*x(1) - 2352.71;
out(3) = 5.67e-8*x(3).^4 + 1.865*x(3) - x(4) - 2250;
out(4) = x(4) - 0.71*x(2) + 7.46*x(3) - 11093;

end

function X = NewtonSystemSol(F,Jac,Xo,tol,maxIt)
xi = Xo(1); yi = Xo(2);
for i = 1:maxIt
   J = Jac(xi,yi);
%    F = -[F1(xi,yi); F2(xi,yi)];
   dp = -J\F(xi,yi);
   xip = xi + dp(1);
   yip = yi + dp(2);
   Err_x = abs((xip - xi)/xi);
   Err_y = abs((yip - yi)/yi);
   
%    fprintf('i = %i  x = %-7.4f  y = %-7.4f  Error in x = %-7.4g Error in y = %-7.4g \n',...
%        i,xip,yip,Err_x,Err_y);
   
   X(1) = xip; X(2) = yip;
   if (Err_x < tol) && (Err_y < tol)
       break
   else
       xi = xip; yi = yip;
   end
    
end
end

function x = Gauss(A,b)
[R,C] = size(A);

% forward elimination phase
for j = 1:R-1 % start with first row, repeat for all but last row
    for i = (j+1):R % start with j+1 row and repeat for all remaining rows
        m = A(i,j)/A(j,j); % calculate pivot
        A(i,j:C) = A(i,j:C) - m*A(j,j:C); % term A(i,j) should now be zero
        % and other members of that row updated accordingly
        b(i) = b(i) - m*b(j); % update the right hand side
    end
end

% back substitution phase
x = nan(R,1);
x(R) = b(R)/A(R,R);
for i = (R-1):-1:1
   x(i) = (b(i) - A(i,(i+1):end)*x((i+1):end))/A(i,i); 
   % The line above is a little tricky. Note the vector addition in the 
   % numerator.
end

% at this point, all values of x should be updated.  That is an assumption
% worth verifying.

% Starting with MATLAB 2022a there is a function called "anynan" to
% determine if any element of an array is equal to nan.  If you have an
% older version of MATLAB, the statement below will also work.

assert(sum(isnan(x))==0,'Error! An element of x is NaN!');

end

function x = GaussPivot(A,b)
[R,C] = size(A);

for j = 1:R-1
   % select pivot from remaining rows in column j
   [~,j_piv] = max(abs(A(j:end,j))); 
   j_piv = j_piv + (j-1); % correct row index returned by max() 
   
   % swap row j and j_piv, b(j) and b(j_piv)
   a_row_tmp = A(j,:); b_tmp = b(j); % save row to be over-written
   A(j,:) = A(j_piv,:); b(j) = b(j_piv);% over-write with new pivot row
   A(j_piv,:) = a_row_tmp; b(j_piv) = b_tmp; % replace row
   
   % introduce zeros in column j below pivot
   for i = (j+1):R
       m = A(i,j)/A(j,j); % calculate pivot
       A(i,j:C) = A(i,j:C) - m*A(j,j:C); % term A(i,j) should now be zero
       % and other members of that row updated accordingly
       b(i) = b(i) - m*b(j); % update the right hand side
   end
end

% back substitution phase (same as with Gauss Elimination without
% Pivoting)
x = nan(R,1);
x(R) = b(R)/A(R,R);
for i = (R-1):-1:1
    x(i) = (b(i) - A(i,(i+1):end)*x((i+1):end))/A(i,i);
    % The line above is a little tricky. Note the vector addition in the
    % numerator.
end

end

