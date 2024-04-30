function [x,y] = BVPShootSecant(fOFx,gOFx,hOFx,a,b,n,Ya,Yb,WL,WH)
% function [x,y] = BVPShootSecant(fOFx,gOFx,hOFx,a,b,n,Ya,Yb,WL,WH)
% solves a second-order BVP of form:
% y'' + f(x)y' + g(x)y = h(x) on the domain a < x < b with 
% y(a) = Ya, and y(b) = Yb.
%
% inputs
% fOFx - function handle
% gOFx - function handle
% hOFx - function handle
% a - scalar domain left boundary
% b - scalar domain right boundary
% n - scalar, number of subdivisions
% Ya - scalar, left BC
% Yb - scalar, right BC
% WL - scalar slope low guess
% WH - scaalr slope high guess for bisection

h = (b-a)/n;

% define two ODEs: x - independent variable, y - vector of 2 dependent
% variables.
ODE1 = @(x,y) y(2);
ODE2 = @(x,y) -fOFx(x).*y(2) - gOFx(x).*y(1) + hOFx(x);

[~,yL,~] = Sys2ODEsRK4(ODE1,ODE2,a,b,h,Ya,WL);
[~,yH,~] = Sys2ODEsRK4(ODE1,ODE2,a,b,h,Ya,WH);

Ei = yL(end) - Yb;
Eim1 = yH(end) - Yb;
Wim1 = WH;
Wi = WL;



%% set bisection iteration parameters
tol = 1e-4;
imax = 100;

for i = 1:imax
   Wip1 = Wi - (Wim1 - Wi)*Ei/(Eim1 - Ei);
   [x,y,~] = Sys2ODEsRK4(ODE1,ODE2,a,b,h,Ya,Wip1);
   Eim1 = Ei;
   Ei = y(end) - Yb;
   
   if abs(Ei) < tol
      break; %<-- success! 
   end
   
   Wim1 = Wi;
   Wi = Wip1;
end

if i == imax
    fprintf('Failed to find a solution within error tolerance in %d iterations\n',...
        imax);
end



