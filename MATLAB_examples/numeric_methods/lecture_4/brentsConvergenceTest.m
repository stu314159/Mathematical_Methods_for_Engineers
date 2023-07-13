%% Brent's Method Convergence Test(approximate)
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
tol = 5e-16;


fprintf('iter       a         b            xNS         f(xNS)         Tol\n');
xNS = brent(F,a,b,maxIt,tol);

function xNS = brent(f,a,b,maxIt,tol)
% ensure that I start with the solution bracketed
fa = f(a); fb = f(b);
assert(a < b, 'Error!  a must be less than b for [a,b]!\n');
assert(fa*fb < 0,'Error!  [a,b] does not bracket a single root!\n');

for j = 1:maxIt
   % do bisection
   c = 0.5*(a + b);
      
   % estimate root with inverse quadratic interpolation (use bisection
   % estimate as point c)
   fc = f(c);
   
   % this is ugly - replace with expression from NR in C.
%    x_iqi = (-fa)*(-fb)*c/((fc-fa)*(fc-fb)) + ...
%        (-fb)*(-fc)*a/((fa - fb)*(fa-fc)) + ...
%        (-fc)*(-fa)*b/((fb - fc)*(fb - fa));

% stil breaks...
   R = fb/fc; S = fb/fa; T = fa/fc;
   P = S*(T*(R-T)*(c-b) - (1-R)*(b-a));
   Q = (T-1)*(R-1)*(S-1);
   x_iqi = b + P/Q;
   
   if (fa*fc < 0)
       a_nb = a;
       b_nb = c; fb = fc;
   elseif (fc*fb < 0)
       a_nb = c; fa = fc;
       b_nb = b;
   else 
       error('The bracket is broken! \n');
   end
   
   % is x_iqi in [a,b]? if not, just go forward with the bisection
   % estimate.
   %if (x_iqi < a) || (x_iqi > b)
   
   % does f(x_iqi) make a smaller bracket than [a,c] or [c,b]?
   % if not, just go forward with the bisection estimate.
   
   % this reduces to: is x_iqi in the right bracket; if it is, then it's
   % better than bisection.
   if (x_iqi > a_nb) && (x_iqi < b_nb)
       fc = f(x_iqi);
       
       if fc == 0 % just in case.
           xNS = x_iqi;
           return;
       end
       
       if(fa*fc < 0)
           a = a_nb;
           b = x_iqi;
       else
           a = x_iqi;
           b = b_nb;
       end    
       
   else
       fprintf('iteration %d used bisection\n',j);
       a = a_nb;
       b = b_nb;
   end
   
   itol = 0.5*(b-a);
   xNS = 0.5*(a+b);
   
   if itol < tol
       return
   end
    
   fprintf('%3i %11.6f %11.6f  %11.6f  %11.6g  %11.6g\n',...
        j,a,b,xNS,f(xNS),itol);
    
end
fprintf('Reached max iterations, best estimate: %11.9g \n',xNS);

end