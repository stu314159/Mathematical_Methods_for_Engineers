%% Example 3.5
% solve using Newton's Method

clear
clc
close 'all'

F1 = @(x,y) y - 0.5*(exp(x./2) + exp(-x./2));
F2 = @(x,y) 9*x.^2 + 25*y.^2 - 225;

dF1x = @(x,y) -0.25*(exp(x./2) - exp(-x./2));
dF1y = @(x,y) 1;

dF2x = @(x,y) 18*x;
dF2y = @(x,y) 50*y;

Jac = @(x,y) [dF1x(x,y) dF1y(x,y); dF2x(x,y) dF2y(x,y)];

xi = 2.5; yi = 2; Err = 1e-10;

for i = 1:5
   J = Jac(xi,yi);
   F = -[F1(xi,yi); F2(xi,yi)];
   dp = J\F;
   xip = xi + dp(1);
   yip = yi + dp(2);
   Err_x = abs((xip - xi)/xi);
   Err_y = abs((yip - yi)/yi);
   
   fprintf('i = %i  x = %-7.4f  y = %-7.4f  Error in x = %-7.4g Error in y = %-7.4g \n',...
       i,xip,yip,Err_x,Err_y);
   
   if (Err_x < Err) && (Err_y < Err)
       break
   else
       xi = xip; yi = yip;
   end
    
end