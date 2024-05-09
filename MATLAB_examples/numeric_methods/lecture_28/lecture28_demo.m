%% Lecture 28 - Solving BVPs with the Shooting Method
clear
clc
close 'all'
%% Define the ODE and all parameters
a = 0; b = 0.1;
Ta = 473; % K, temperature at base of fin
Tb = 293; % K, temperature at the tip of the fin
N = 1000; %
ode = @(x,T) pin_fin(x,T);
x = linspace(a,b,N);

%% set RK4 solver tableau parameters
% note: any other ODE solver would also suffice.
solver = @odesExplicitRK;
s = 4;
BT = zeros(s+1,s+1); % Butcher Tableau for RK method
C = [0; 1/2; 1/2; 1; 0]; % sample points
B = [0 1/6 2/6 2/6 1/6]; % weights
A = [0 0 0 0;    % RK matrix
    1/2 0 0 0;
    0 1/2 0 0;
    0 0 1 0;];
BT(:,1) = C;
BT(end,:) = B;
BT(1:s,2:end) = A;

%% Commence Shooting Method
%  Function that will take the ODE, a solver, and a guessed value for the
%  shooting method
odeWithSolver = @(dTa_g) solver(ode,a,b,N,[Ta;dTa_g],BT);

%  Function that will solve the ODE with guessed value and return 0 when
%  the "correct" value is guessed.  This can be passed to any root-finding
%  function.
tgt_err_fun = @(dTa_g) shot_function(odeWithSolver,dTa_g) - Tb;

% need two guesses at the slope
dT_a_H = 0; % slope = 0 at x=a means no heat transfer out.
dT_a_L = 4*(Tb - Ta)/(b-a); % guess a strong linear function

method = 4;
tol = 1e-4; % tolerance on convergence to the BC
switch method
    case 1
        dT_a = BisectionRoot(tgt_err_fun,dT_a_H,dT_a_L,tol);
    case 2
        imax = 100;
        dT_a = SecantRoot(tgt_err_fun,dT_a_H,dT_a_L,tol,imax);
        
    case 3
        dT_a = fzero(tgt_err_fun,dT_a_H);
        
    case 4
        dT_a = fsolve(tgt_err_fun,dT_a_H);
        
    otherwise
        error('Invalid Solver Choice!');
end
%% Solve one last time with converged dT_a
T = odeWithSolver(dT_a);

figure(1)
plot(x,T(1,:),'-b','linewidth',3);
title('Pin Fin Temperature Profile',...
    'fontsize',14,'fontweight','bold');
xlabel('X [m]','fontsize',12,'fontweight','bold');
ylabel('T [^\circ C]','fontsize',12,'fontweight','bold');
grid on
set(gca,'fontsize',10,'fontweight','bold');
%% Local Functions
function Tb = shot_function(odeWithSolver,dT_a)
% given an estimate for T'(a), solves the IVP and returns
% the value of T(b).
T_trial = odeWithSolver(dT_a);
Tb = T_trial(1,end);
end

function F = pin_fin(~,T)
% define the IVP
% two args so it can work with ode45
Ac = 1.6e-5; % m^2, fin cross sectional area
P = 0.016; % m, perimeter of pin cross section
h_c = 40; % W/m^2-K, convective heat transfer coefficient of air around pin
k = 250; % W/m-K, thermal conductivity of pin material
emiss = 0.5; % emissivity of pin material
sigma_sb = 5.67e-8; % W/m^2-K^4, Stefan-Boltzmann constant
Ts = 293; % K, temperature of surrounding air

F = nan(2,1);
F(1) = T(2);
F(2) = ((h_c*P)/(k*Ac))*(T(1) - Ts) + ...
    ((emiss*sigma_sb*P)/(k*Ac))*(T(1).^4 - Ts^4);

end

function y = odesExplicitRK(ODE,a,b,N,yINI,BT)
% function y = odeExplicitRK(ODE,a,b,h,yINI,BT)
% y = solution (vector)
% ODE = function handle for y'
% a,b = begining and end of the interval for solution
% N = number of steps between a and b
% yINI = initial value for the solution
% BT = Butcher Tableau

% get Butcher Tableau Parameters
s = length(BT)-1;
c = BT(1:s,1);
B = BT(s+1,2:end);
A = BT(1:s,2:end);
stages = s;

x = linspace(a,b,N);
sys_size = length(yINI);
y = nan(sys_size,N);
y(:,1) = yINI;
h = x(2)-x(1);
for t = 1:(N-1)
    Xi = nan(sys_size,stages);
    
    for s = 1:stages
       Xi(:,s) = y(:,t);
       for i = 1:(s-1)
          Xi(:,s) = Xi(:,s) + h*A(s,i)*ODE(x(t)+c(i)*h,Xi(:,i)); 
       end
    end
    
    y(:,t+1) = y(:,t);
    for i = 1:stages
       y(:,t+1) = y(:,t+1) + h*B(i)*ODE(x(t)+c(i)*h,Xi(:,i)); 
    end
    
end
end

function x_bi = BisectionRoot(fun,a,b,TOL)

tol = TOL;
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

function Xs = SecantRoot(Fun,Xa,Xb,Err,imax)

for i = 1:imax
   FXb = Fun(Xb);
   Xi = Xb - FXb*(Xa - Xb)/(Fun(Xa) - FXb);
   
   if abs((Xi - Xb)/Xb) < Err
       Xs = Xi;
       break;
   end
   Xa = Xb;
   Xb = Xi;
     
        
end

end