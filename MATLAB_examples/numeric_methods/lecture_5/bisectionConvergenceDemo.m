%% Bisection Convergence Demo

clear
clc
close 'all'

% find the square root of two
F = @(x) x.^2 - 2;

a = 0; b = 2; % initial bracket
Fa = F(a); Fb = F(b);
% verify that Fa*Fb > 0
if (Fa*Fb > 0)
   error('Error: The function has the same sign at points a and b'); 
end

Err = 5e-16;

imax = 150;
fprintf('Exact solution to full double precetion = %16.15f \n',sqrt(2));


fprintf('iter       a         b            xNS         f(xNS)         Tol\n');
for i = 1:imax
    xNS = (a + b)/2;
    toli = (b - a)/2;
    FxNS = F(xNS);
    fprintf('%3i %11.6f %11.6f  %11.6f  %11.6f  %11.6g\n',...
        i,a,b,xNS,FxNS,toli);
    
    if FxNS == 0
        fprintf('An exact solution x = %11.6f was found\n',xNS);
        break;
    end
    
    if abs(toli/sqrt(2)) < Err % check relative error
        fprintf('Success!! Xi = %16.15f \n',xNS);
        break; % end the loop since I meet my error tolerance
    end
    
    % if toli > tol and we've reached maximum iterations, give up.
    if i == imax
        fprintf('Solution was not obtained in %i iterations\n',imax);
        break
    end
    
    % update bracket
    if (F(a)*FxNS < 0)
        b = xNS;
    else
        a = xNS;
    end
    
end