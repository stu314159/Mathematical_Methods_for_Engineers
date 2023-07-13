%% Lecture 3 - Bisection Method Example
% Implementation roughly based on example 3-1 of  your textbook

clear
clc
close 'all'

e = 0.0015; % mm
D = 4.0; % mm
Re = 13700;

F = @(f) 1./sqrt(f) + 2.0*log10((e/D)/3.7 + 2.51/(Re*sqrt(f)));

a = 0.001; b = 1; 
imax = 20; tol = 0.0001;

Fa = F(a); Fb = F(b);

% verify that Fa*Fb > 0

if (Fa*Fb > 0)
   error('Error: The function has the same sign at points a and b'); 
end

fprintf('iter       a         b            xNS         f(xNS)         Tol\n');

for i = 1:imax
   xNS = (a + b)/2;
   toli = (b - a)/2;
   FxNS = F(xNS);
   fprintf('%3i %11.6f %11.6f  %11.6f  %11.6f  %11.6f\n',...
       i,a,b,xNS,FxNS,toli);
   
   if FxNS == 0
       fprintf('An exact solution x = %11.6f was found\n',xNS);
       break;
   end
   
   if toli < tol
       fprintf('Success!! x = %11.6f \n',xNS);
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