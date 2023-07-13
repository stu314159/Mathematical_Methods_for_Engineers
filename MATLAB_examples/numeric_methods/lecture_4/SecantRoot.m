function xNS = SecantRoot(F,Xa,Xb,Err,imax)
% SecantRoot finds the root of F=0 near the points Xa and Xb using the Secant
% Method
% Input Variables:
% F       Name of a user-defined function that calculates F(x) for a given x
% Xa, Xb  Two points in the vicinity of the desired root
% Err     Estimated Relative Error at which iterations should be terminated
% imax    maximum number of iterations
%
% Output Variables:
% xNS     Solution

for i = 1:imax
    FXb = F(Xb);
    Xi = Xb - FXb*(Xa - Xb)/(F(Xa) - FXb);
    
    if abs((Xi - Xb)/Xb) < Err
        xNS = Xi;
        break;
    end
    
    Xa = Xb;
    Xb = Xi;
end

if i == imax
    fprintf('Solution was not obtained in %i iterations. \n',imax);
    xNS = nan;
end