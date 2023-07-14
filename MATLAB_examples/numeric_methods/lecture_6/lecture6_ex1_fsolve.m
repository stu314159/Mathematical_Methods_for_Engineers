%% Use fsolve for system of nonlinear equations
clear
clc
close 'all'

F = @(x) ex3p5(x);
X0 = [2.5, 2.0];

option_set = 1;
% 1 = no outputs
% 2 = detailed outputs

switch option_set
    
    case 1
        options = optimoptions('fsolve','Display','none');
        
    case 2
        options = optimoptions('fsolve','Display','iter-detailed',...
            'MaxIterations',1000,'StepTolerance',1e-10);
        
    otherwise % default options
        options = optimoptions('fsolve');
        
end

[x,fval,exitflag,output] = fsolve(F,X0,options);

fprintf('Root found at x = %8.7g, y = %8.7g \n',x(1),x(2));
fprintf('fval = \n'); disp(fval);
fprintf('exitflag = %d \n',exitflag);

%% Local Function
function out = ex3p5(x)
[m,n] = size(x); % expect scalar or vector input
assert(min(m,n)==1,'Error!  Vector input expected for x! \n');
out = nan(m,n); % construct output

out(1) = x(2) - 0.5*(exp(x(1)./2) + exp(-x(1)./2));
out(2) = 9*x(1).^2 + 25*x(2).^2 - 225;

end