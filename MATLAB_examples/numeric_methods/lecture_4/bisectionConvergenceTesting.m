%% Bisection Convergence Testing

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

Fa = F(a); Fb = F(b);
% verify that Fa*Fb > 0
if (Fa*Fb > 0)
   error('Error: The function has the same sign at points a and b'); 
end

Err = 2.2204e-16;

imax = 150;



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