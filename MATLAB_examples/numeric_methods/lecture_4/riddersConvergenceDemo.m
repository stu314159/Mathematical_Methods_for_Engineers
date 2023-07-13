%% Ridders' Method Convergence Demo
% From NR in C section 9.2 pg 358

% it seems to work.  I confess that I really do not understand the
% implementation.  Transcribed from Numerical Recipes.

clear
clc
close 'all'

% find the square root of two
F = @(x) x.^2 - 2;
x1 = 0; x2 = 2; % initial bracket

Err = 1e-16;
imax = 150;

fprintf('Exact solution to full double precision = %16.15f \n',sqrt(2));
x = sqrt(2);

fprintf('iter       a         b            xNS         f(xNS)         Tol\n');
xNS = ridder(F,x1,x2,imax,Err);

fprintf('Root found: %16.15f \n',xNS);


function xNS = ridder(F, x1, x2, MAXIT, xacc)
fL = F(x1); fH = F(x2);

% verify that x1 and x2 bracket a root
assert (fL*fH < 0,"[x1,x2] does not bracket a single root!\n");
xL = x1; xH = x2;

% check the crazy case where one or the other end point is a root
if fL == 0
    xNS = x1;
    return
elseif fH == 0
    xNS = x2;
    return
end

ANS = nan; % initialize the return value

for j = 1:MAXIT
    
   xm = 0.5*(xL + xH);
   fm = F(xm);
   s = sqrt(fm*fm - fL*fH);
   if s == 0
       xNS = ANS;
       return
   end
   
   if fL >= fH
       fact1 = 1;
   else 
       fact1 = -1;
   end
   xnew = xm+(xm-xL)*fact1*fm/s;
   
   if abs(xnew - ANS) <= xacc
       xNS = ANS;
       return;
   end
   ANS = xnew;
   fnew = F(ANS);
   
   if fnew == 0
       xNS = ANS;
       return;
   end
   
   % bookkeeping to keep the root bracketed on next iteration
   if (SIGN(fm,fnew) ~= fm)
       xL = xm;
       fL = fm;
       xH = ANS;
       fH = fnew;
   elseif(SIGN(fL,fnew) ~= fL)
       xH = ANS;
       fH = fnew;
   elseif(SIGN(fH,fnew) ~= fH)
       xL = ANS;
       fL = fnew;
   else
       error('Something went wrong!\n');
   end
   
   toli = abs(xH - xL);
   fprintf('%3i %11.6f %11.6f  %11.6f  %11.6g  %11.6g\n',...
        j,xL,xH,ANS,fnew,toli);
   
   
   if(toli <= xacc)
       xNS = ANS;
       return;
   end
   
   
   
end
error('Exceeded maximum iterations!\n');

end

function y = SIGN(a,b)
y = abs(a)*sign(b);
end