%% fzero Convergence Test(approximate)
clear
clc
close 'all'

func = 2;
switch func  
    case 1
        F = @(x) sqrt(x) - cos(x);
        a = 0; b = 1;
        
    case 2
        F = @(x) exp(x) - x.^2 + 3*x - 2;
        a = 0; b = 1;
        
    case 3
        F = @(x) x.^3 - 7*x.^2 + 14*x - 6;
        a = 3.2; b = 4.0;
    
end

maxIt = 150;
tol = 2.2204e-16;


options = optimset('Display','iter','TolX',tol);
[x,fval,exitflag,output] = fzero(F,[a,b],options);
fprintf('root: %16.15f \n',x);
fprintf('final fval: %16.15g \n',fval);
