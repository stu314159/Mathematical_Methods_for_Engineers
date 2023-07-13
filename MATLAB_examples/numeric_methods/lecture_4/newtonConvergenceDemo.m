%% Newton Convergence Demo

clear
clc
close 'all'

% find the square root of two
F = @(x) x.^2 - 2;
dF = @(x) 2*x;

Xest = 1.; % initial guess

Err = 1e-16;
imax = 5;
fprintf('Exact solution to full double precetion = %16.15f \n',sqrt(2));
for i = 1:imax
   Xi = Xest - F(Xest)/dF(Xest); 
    
   if abs((Xi - Xest)/Xest) < Err
       xNS = Xi;
       break;
   end
   
   fprintf('After %i iterations, Xi = %16.15f \n',i,Xi);
   
   Xest = Xi;
      
end