
clear
clc
close all


% Fix sample point xI(2) at 0.5. (on domain [0,1]) 
% Mathematica used to find: xI(1),xI(3) and weight points wI(1)-wI(3)
% to exactly solve any 4th order polynomial.

% As seen by the associated Mathematic Notebook, if one were to free the
% xI(2) and look for a 3-point 5th order integration scheme, this is what
% you would get. Interesting.

% For this implementation, the relative error is higher than it should be
% due to the limited precision of the sample points and weights.
xI = NaN(3,1);
xI(1) = 0.1127016653792583;
xI(2) = 0.5;
xI(3) = 0.887298334620755;

wI = NaN(3,1);
wI(1) = 0.277777777777777786;
wI(2) = 0.444444444444444493;
wI(3) = wI(1);

f = @(x) exp(-(x.^2));
a = -3;
b = 3;

% refine so that the integration scheme handles arbitrary 1-D domains.
Jac = (b-a);
xTrans = @(t) a + (b-a)*t;

intExact = integral(f,a,b);

intApprox1 = f (xTrans(xI))'*wI*Jac;

RelErr = norm (intExact - intApprox1,2)/intExact;
fprintf('Relative Error with 1 subinterval = %g.\n',RelErr);

% make into a composite rule:

N = 40;

% Split the domain into N subintervals; execute the integration scheme on
% each of the sub-intervals and add the result.
xSplit = linspace(a,b,N+1); 
Jac = xSplit(2)-xSplit(1);
% not very efficient, but for this test implementation, simply loop through
% all of the subdomains, take the integration and sum the result:

intApproxN = 0;
for d = 1:N
   xMin = xSplit(d);
   xMax = xSplit(d+1);
   xT = xMin + (xMax-xMin)*xI;
   intApproxN = intApproxN + f(xT)'*wI*Jac;
end

RelErrN = norm(intExact - intApproxN,2)/intExact;
fprintf('Relative Error with %i subintervals = %g.\n',N,RelErrN);

%% Convergence Rate of 3-point 6th order method
Nspace = 2.^(3:7);
intF_vec = NaN(1,length(Nspace));
relError_vec = intF_vec;
for n = 1:length(Nspace)
    intF_vec(n) = myThreePointFifthOrder(f,a,b,Nspace(n));
    relError_vec(n) = abs(intF_vec(n) - intExact)/intExact;
end
loglog(Nspace,relError_vec,'-b','linewidth',2);
grid on
hold on
loglog(Nspace,.1*(1./Nspace).^6,'-.r','linewidth',2)
title('Convergence of 3-point 6th Order Rule','fontsize',16,'fontweight','bold');
xlabel('Number of Subdivisions','fontsize',14,'fontweight','bold');
ylabel('Relative Error','fontsize',14,'fontweight','bold');
legend('Relative Error','Sixth Order Convergence');
set(gca,'fontweight','bold');


