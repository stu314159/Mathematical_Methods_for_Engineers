function xNS = NewtonRoot(F,dF,Xest,Err,imax)
% Newton Root finds the root of F=0 near the point Xest using Newton's
% Method
% Input Variables:
% F     Name of a user-defined function that calculates F(x) for a given x
% dF    Name of a user-defined function that calculates F'(x) for a given x
% Xest  Initial estimate of the solution
% Err   Estimated Relative Error at which iterations should be terminated
% imax  maximum number of iterations
%
% Output Variable:
% xNS   solution

for i = 1:imax
   Xi = Xest - F(Xest)/dF(Xest);
   if abs((Xi - Xest)/Xest) < Err
       xNS = Xi;
       break;
   end
   
   Xest = Xi;
    
end

if i == imax
    fprintf('Solution was not obtained in %i iterations. \n',imax);
    xNS = nan;
end