%% BiRegula Falsi Convergence Test
% this is a prototype to a solution for problem 3.21 in the text.

clear
clc
close 'all'

func = 1;
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
xNSm = b; % just to get started...
% verify that Fa*Fb > 0
if (Fa*Fb > 0)
   error('Error: The function has the same sign at points a and b'); 
end

Err = 2.2204e-16;

imax = 150;


fprintf('iter       a         b            xNS         f(xNS)         ERE\n');
for i = 1:imax
    % bisection first:
    xNSb = (a + b)/2;
    FxNS = F(xNSb);
    toli = (b - a)/2;
    % update bracket
    if (F(a)*FxNS < 0)
        b = xNSb;
    else
        a = xNSb;
    end
    
    % now do RF 
    xNS = (a*F(b) - b*F(a))./(F(b) - F(a));
   
    FxNS = F(xNS);
    
    ERE = abs((xNS - xNSm)/xNSm); % estimated relative error
    
    fprintf('%3i %11.6f %11.6f  %11.6f  %11.6f  %11.6g\n',...
        i,a,b,xNS,FxNS,ERE);
    
%     if FxNS == 0
%         fprintf('An exact solution x = %16.15f was found\n',xNS);
%         break;
%     end
    
    if ERE < Err % Use True Relative Error
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
    xNSm = xNS; % save this for the next iteration.
end