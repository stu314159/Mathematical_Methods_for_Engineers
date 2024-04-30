%% Assignment #1 "Solution" script
clear
clc
close 'all'

%% Problem 3 (1.4)
fprintf('\n\n\n Problem #3 \n\n\n');
n = 10;
[a_pi_a,tre_a] = approxPi(n);
fprintf('N = %d: Approximate value: %g, True Relative Error: %g \n',...
    n,a_pi_a,tre_a);

n = 20;
[a_pi_b,tre_b] = approxPi(n);
fprintf('N = %d: Approximate value: %g, True Relative Error: %g \n',...
    n,a_pi_b,tre_b);

n = 40;
[a_pi_c,tre_c] = approxPi(n);
fprintf('N = %d: Approximate value: %g, True Relative Error: %g \n',...
    n,a_pi_c,tre_c);

%% Problem 5
fprintf('\n\n\n Problem #5 \n\n\n');
f = @(x) x - 2*exp(-x);
a = 0; b = 1;

x_bi = myBisection(f,a,b);
fprintf('Estimated root: %g \n',x_bi);

%% Problem 6
fprintf('\n\n\n Problem #6 \n\n\n');
R = 0.08206; %(L-atm/mole-K), gas constant
a = 1.39; % L^2-atm/mole^2, constant
b = 0.03913; % L/mole, constant
n = 1.5; % moles
T = 25 + 273; % K, temperature of the gas.
P = 13.5; % atm, pressure

f = @(V) (n*R*T)./(V-n*b) - (n^2*a)./(V.^2) - P;
%% Just to get an idea where the zeros are.
fplot(f,[1 100])
%%
a_bi = 1; b_bi = 10; % could plot first to get an idea 
% where the zero lies.
V_bi = myBisection(f,a_bi,b_bi);

fprintf('Volume = %g liters. \n',V_bi);


%% Local Functions
function [ap_pi, tre] = approxPi(n)
ap_pi = 0;
for k = 1:n
    ap_pi = ap_pi + 4*((-1)^(k-1))*(1/(2*k-1));
end

tre = abs(pi - ap_pi)/pi;

end

function x_bi = myBisection(fun,a,b)

tol = 1e-6;
x_bi = (a+b)/2;
% verify that fun(a) and fun(b) have different sign
if(fun(a)*fun(b)>0)
   error('No root lies between %g and %g!\n',a,b); 
end

maxIt = 1e2;
for k = 1:maxIt
    FxNS = fun(x_bi);
    
    % check if FxNS is within the tolerance
    if (abs(FxNS) < tol)
        return;
    end
    
    % update brackets
    if (fun(a)*FxNS < 0)
        b = x_bi;
    else
        a = x_bi;
    end
    
    % update estimated x_bi based on new bracket
    x_bi = (a+b)/2;
       
end

% if I ever get here, then I failed to meet the tolerance 
% within the maximum number of iterations
fprintf('Warning: failed to find root within specified tolerance.\n');
fprintf('Current root: %g, fun(x_bi) = %g. \n',x_bi,fun(x_bi));


end