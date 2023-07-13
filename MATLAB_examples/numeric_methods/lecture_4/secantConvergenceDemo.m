%% Secant Convergence Demo

clear
clc
close 'all'

% find the square root of two
F = @(x) x.^2 - 2;

Xa = 0; 
Xb = 2;

Err = 1e-16;
imax = 15;
fprintf('Exact solution to full double precetion = %16.15f \n',sqrt(2));
for i = 1:imax
   FXb = F(Xb);
   Xi = Xb - FXb*(Xa - Xb)/(F(Xa) - FXb);
   
   if abs((Xi - Xb)/Xb) < Err
       xNS = Xi;
       break;
   end
   Xa = Xb;
   Xb = Xi;
   
   fprintf('After %i iterations, Xi = %16.15f \n',i,Xi);
   
        
end