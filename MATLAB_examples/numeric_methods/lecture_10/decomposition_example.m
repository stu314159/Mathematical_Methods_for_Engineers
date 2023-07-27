clear
clc
close 'all'

n = 5;
A = rand(n,n);
b1 = rand(n,1); b2 = rand(n,1);

Ad = decomposition(A,'auto');

x1 = Ad\b1;
x2 = Ad\b2;


rel_resid = norm(A*x1-b1,2)/norm(b1,2);
fprintf('Relative residual = %g \n',...
    rel_resid);