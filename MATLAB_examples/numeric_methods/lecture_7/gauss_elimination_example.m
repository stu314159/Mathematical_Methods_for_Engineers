%% Example: Gauss Elimination
clear
clc
close 'all'

%% System to be solved
A = [4, -2, -3, 6;
    -6, 7, 6.5, -6;
    1, 7.5, 6.25, 5.5;
    -12, 22, 15.5, -1];

b = [12;
    -6.5;
    16;
    17];

x = Gauss(A,b);

fprintf('Solution with local function Gauss: \n');
disp(x);

fprintf('Solution with backslash: \n');
x_gold = A\b;
disp(x_gold);

rel_error = norm(x-x_gold,2)/norm(x_gold,2);
fprintf('Relative difference: %g \n',rel_error);

%% Local function implementing Gauss Elimination
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
   % The line above is a little tricky. Note the vector multiplication 
   % in the numerator.
end

% at this point, all values of x should be updated.  That is an assumption
% worth verifying.

% Starting with MATLAB 2022a there is a function called "anynan" to
% determine if any element of an array is equal to nan.  If you have an
% older version of MATLAB, the statement below will also work.

assert(sum(isnan(x))==0,'Error! An element of x is NaN!');

end
