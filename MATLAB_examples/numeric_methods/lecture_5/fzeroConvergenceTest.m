%% fzero Convergence Test(approximate)
clear
clc
close 'all'

func = 1;
switch func  
    case 1
        F = @(x) x.^2 - 2;
        a = 0; b = 2;
        
    case 2
        F = @(x) exp(x) - x.^2 + 3*x - 2;
        a = 0; b = 1;
        
    case 3
        F = @(x) x.^3 - 7*x.^2 + 14*x - 6;
        a = 3.2; b = 4.0;
    
end

maxIt = 150;
tol = 5e-16;


options = optimset('Display','iter','TolX',tol);
%options = optimset('TolX',tol,'PlotFcns',@optimplotx);
%options = optimset('TolX',tol','PlotFcns',@optimplotfval);
[x,fval,exitflag,output] = fzero(F,[a,b],options);
fprintf('root: %16.15f \n',x);
fprintf('final fval: %16.15g \n',fval);
