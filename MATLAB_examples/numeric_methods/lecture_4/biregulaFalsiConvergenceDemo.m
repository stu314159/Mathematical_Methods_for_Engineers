%% BiRegula Falsi Convergence Demo
% this is a prototype to a solution for problem 3.21 in the text.

clear
clc
close 'all'

% find the square root of two
F = @(x) x.^2 - 2;

a = 0; b = 2; % initial bracket
xNSm = b; % just to get us started.

Fa = F(a); Fb = F(b);
% verify that Fa*Fb > 0
if (Fa*Fb > 0)
   error('Error: The function has the same sign at points a and b'); 
end

Err = 1e-16;

imax = 150;
fprintf('Exact solution to full double precetion = %16.15f \n',sqrt(2));
x = sqrt(2);

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
    TRE = abs(xNS - x)/abs(x);
    ERE = abs((xNS - xNSm)/xNSm); % estimated relative error
    
    fprintf('%3i %11.6f %11.6f  %11.6f  %11.6f  %11.6g\n',...
        i,a,b,xNS,FxNS,ERE);
    
    if FxNS == 0
        fprintf('An exact solution x = %11.6f was found\n',xNS);
        break;
    end
    
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