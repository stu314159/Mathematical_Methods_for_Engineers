%% Lecture 27 - MATLAB Demo
clear
clc
close 'all'

%% Problem to be solved
F = @(t,y) (t.^2).*y;
yINI = 1;

xMin = 0; xMax = 1;
%% Exact solution
y_gold = @(t) exp((1/3)*(t.^3));
x_gold = linspace(xMin,xMax,1000);

figure(1)
plot(x_gold,y_gold(x_gold),'-c','linewidth',3);
grid on
title('Exact Solution','fontsize',14,'fontweight','bold');
xlabel('X','fontsize',12,'fontweight','bold');
ylabel('Y(X)','fontsize',12,'fontweight','bold');

%% Solve with ODE23
myRT = 1e-10;
options = odeset('RelTol',myRT);
[t_ode23,y_ode23] = ode23(F,[xMin xMax],yINI,options);

figure(2)
plot(t_ode23,y_ode23,'-sr','linewidth',3,'markersize',4);
title('Solution with ODE23','fontsize',16,'fontweight','bold');
grid on
xlabel('T','fontsize',14,'fontweight','bold');
ylabel('Y(T)','fontsize',14,'fontweight','bold');
set(gca,'fontsize',12,'fontweight','bold');

%% Solve with Bogacki-Shampine Embedded 2/3 Method
c = [0 1/2 3/4 1 2 3]'; % sample points;
% encode the order of embedded methods in first column, last two rows.
A = [0 0 0 0;
    1/2 0 0 0;
    0 3/4 0 0;
    2/9 1/3 4/9 0]; % RK matrix
b3 = [2/9 1/3 4/9 0]; % weights for 3rd order method
b2 = [7/24 1/4 1/3 1/8]; % weights for 2nd order method
EBT = zeros(6,5);
EBT(:,1) = c;
EBT(1:4,2:end)=A;
EBT(5,2:end) = b2;
EBT(6,2:end) = b3;
[t_bs23,y_bs23] = embeddedRK(F,[xMin xMax],yINI,myRT,EBT);

% get relative error measure at solved points
rel_err_bs23 = norm(y_gold(t_bs23)-y_bs23,2)/...
    norm(y_gold(t_bs23),2);
fprintf('Relative error BS23: %g \n',rel_err_bs23);

%% Plot the result
figure(3)
plot(t_bs23,y_bs23,'-sr','linewidth',3,'markersize',4);

title('Solution with BS23','fontsize',16,'fontweight','bold');
grid on
xlabel('T','fontsize',14,'fontweight','bold');
ylabel('Y(T)','fontsize',14,'fontweight','bold');
set(gca,'fontsize',12,'fontweight','bold');

%% Solve with Fehlberg 4/5
c = [0 1/4 3/8 12/13 1 1/2 4 5]'; % sample points;
% encode the order of embedded methods in first column, last two rows.
A = [0 0 0 0 0 0;
    1/4 0 0 0 0 0;
    3/32 9/32 0 0 0 0;
    1932/2197 -7200/2197 7296/2197 0 0 0;
    439/216 -8 3680/513 -845/4104 0 0;
    -8/27 2 -3544/2565 1859/4104 -11/40 0]; % RK matrix
b5 = [16/135 0 6656/12825 28561/56430 -9/50 2/55]; % weights for 5nd order method
b4 = [25/216 0 1408/2565 2197/4104 -1/5 0]; % weights for 4rd order method
EBT = zeros(8,7);
EBT(:,1) = c;
EBT(1:6,2:end)=A;
EBT(7,2:end) = b4;
EBT(8,2:end) = b5;

[t_rkf45,y_rkf45] = embeddedRK(F,[xMin xMax],yINI,myRT,EBT);

rel_err_rkf45 = norm(y_gold(t_rkf45)-y_rkf45,2)/...
    norm(y_gold(t_rkf45),2);
fprintf('Relative error RKF45: %g \n',rel_err_rkf45);

%% Plot the result
figure(4)
plot(t_rkf45,y_rkf45,'-sr','linewidth',3,'markersize',4);
title('Solution with RKF45','fontsize',16,'fontweight','bold');
grid on
xlabel('T','fontsize',14,'fontweight','bold');
ylabel('Y(T)','fontsize',14,'fontweight','bold');
set(gca,'fontsize',12,'fontweight','bold');

%% Solve with ODE45
options = odeset('RelTol',myRT);
[t_ode45,y_ode45] = ode45(F,[xMin xMax],yINI,options);


%% Plot the result
figure(5)
plot(t_ode45,y_ode45,'-sr','linewidth',3,'markersize',4);
title('Solution with ODE45','fontsize',16,'fontweight','bold');
grid on
xlabel('T','fontsize',14,'fontweight','bold');
ylabel('Y(T)','fontsize',14,'fontweight','bold');
set(gca,'fontsize',12,'fontweight','bold');

%% Solve with ODE45
% re-run without additional points for a "beautiful plot"
options = odeset('RelTol',myRT,'Refine',1,'MaxStep',1);
[t_ode45,y_ode45] = ode45(F,[xMin xMax],yINI,options);

rel_err_ode45 = norm(y_gold(t_ode45)-y_ode45,2)/...
    norm(y_gold(t_ode45),2);
fprintf('Relative error ODE45: %g \n',rel_err_ode45);

%% Plot the result
figure(6)
plot(t_ode45,y_ode45,'-sr','linewidth',3,'markersize',4);
title('Solution with ODE45','fontsize',16,'fontweight','bold');
grid on
xlabel('T','fontsize',14,'fontweight','bold');
ylabel('Y(T)','fontsize',14,'fontweight','bold');
set(gca,'fontsize',12,'fontweight','bold');

%% Solve higher order eqn with ODE45
F = @(t,y) ex2(t,y);
yINI = [-1 2]; % initial values
xMin = 0; xMax = 2.1;
y_gold = @(x) exp(-x./2).*(-cos(2*x)+0.75*sin(2*x));

options = odeset('RelTol',myRT,'Refine',1,'MaxStep',1);
[t_ode45,y_ode45] = ode45(F,[xMin xMax],yINI,options);

rel_err_ode45 = norm(y_gold(t_ode45)-y_ode45(:,1),2)/...
    norm(y_gold(t_ode45),2);
fprintf('Relative error ODE45: %g \n',rel_err_ode45);
% solutions are column vectors
%% Plot the result
figure(7)
plot(t_ode45,y_ode45(:,1),'-sr','linewidth',3,'markersize',4);
title('Solution with ODE45','fontsize',16,'fontweight','bold');
grid on
xlabel('T','fontsize',14,'fontweight','bold');
ylabel('Y(T)','fontsize',14,'fontweight','bold');
set(gca,'fontsize',12,'fontweight','bold');

%% Try embedded RK with new problem.
[t_rkf45,y_rkf45] = embeddedRK(F,[xMin xMax],yINI,myRT,EBT);

rel_err_rkf45 = norm(y_gold(t_rkf45)-y_rkf45(1,:),2)/...
    norm(y_gold(t_rkf45),2);
fprintf('Relative error RKF45: %g \n',rel_err_rkf45);

%% Plot the result
figure(8)
plot(t_rkf45,y_rkf45(1,:),'-sr','linewidth',3,'markersize',4);
title('Solution with RKF45','fontsize',16,'fontweight','bold');
grid on
xlabel('T','fontsize',14,'fontweight','bold');
ylabel('Y(T)','fontsize',14,'fontweight','bold');
set(gca,'fontsize',12,'fontweight','bold');

%% Use MATLAB ODE object (R2023b and later)

M = ode;
M.ODEFcn = F;
M.InitialValue = yINI;
M.RelativeTolerance = myRT; 

sol_M_ode = solve(M,xMin,xMax);

rel_err_ODE = norm(y_gold(sol_M_ode.Time)-sol_M_ode.Solution(1,:),2)/...
    norm(y_gold(sol_M_ode.Time),2);

fprintf('Relative error ODE: %g \n',rel_err_ODE);
%% Plot the result

figure(9)
plot(sol_M_ode.Time,sol_M_ode.Solution(1,:),...
    '-sr','linewidth',3,'markersize',4);
title('Solution with ODE','fontsize',16,'fontweight','bold');
grid on
xlabel('T','fontsize',14,'fontweight','bold');
ylabel('Y(T)','fontsize',14,'fontweight','bold');
set(gca,'fontsize',12,'fontweight','bold');


%% Local Functions
function [t,y] = embeddedRK(F,tspan,y0,RTOL,EBT)
t_sz = 5000; % initial size for output arrays
tsMax = 100000;
sys_size = length(y0);
t = nan(1,t_sz);
y = nan(sys_size,t_sz);

% set initial values
t(1) = tspan(1);
y(:,1) = y0;

[m,~] = size(EBT);
s = m-2;
C = EBT(1:s,1);
A = EBT(1:(end-2),2:end);
BW = EBT(s+1,2:end);
BZ = EBT(s+2,2:end);
p = EBT((end-1),1);

SF = 0.99; % "safety factor"
theta = 1e-14; % factor protects against small values of w.

tStart = tspan(1);
tEnd = tspan(2);
h = (tEnd-tStart)/10; % initial step size
h_new = h;
stopFlag = 0;

for ts = 1:tsMax
    cT = t(ts); % current time
    % find acceptable time step size
    int_it_count = 0; % limit iterations in error control.
    while 1
        int_it_count = int_it_count + 1;
        h = h_new; % update with new step-size
        
        % if cT + h > tEnd, reduce h so we stop right on time
        if cT+h > tEnd
           h = tEnd-cT; 
           stopFlag = 1; % stop after this time step
        end
        
        y(:,ts+1) = y(:,ts);
        K = getSlopeEst(F,t(ts),y(:,ts),h,C,A);
        [w,z] = getWZ(y(:,ts),h,K,BW,BZ);
        err_ts = abs(w-z); % error vector (length = # dofs)
        rel_err_ts = err_ts./max(abs(w),theta); % relative error vector
        max_rel_err_ts = norm(rel_err_ts,inf);
        if max_rel_err_ts < RTOL % time step accepted, compute new (bigger) h
            h_new = SF*(RTOL/max_rel_err_ts)^(1/(p+1))*h;  
            y(:,ts+1) = z; % local extrapolation (use best solution)
            t(ts+1) = cT + h; % update current time   
            break; % exit the while loop
        else
            h_new = h/2;
            if stopFlag == 1
                % unset the stop flag if it was previously set with
                % larger time step.
                stopFlag = 0;
            end
        end % if max_rel_err_ts...
        
        if int_it_count > 10
            error('Error-control is broken!'); % stop the program.
        end
        
        
    end % while 1
    
    % here I should have updated solution and step size for next time step.
    if stopFlag == 1 % end of time stepping.  exit loop.
        break;
    end
end %ts
% trim the output variables.
y = y(:,1:(ts+1));
t = t(1:(ts+1));
end % function

function K = getSlopeEst(F,t,y,h,C,A)
sys_size=length(y);% gen # of dep vars
[s,~] = size(A); % get number of stages
K = zeros(sys_size,s);
K(:,1) = F(t,y);
for i = 2:s
    y_it = y(:);
    x_it = t + C(i)*h;
    for j = 1:(i-1)
        y_it = y_it + h*A(i,j)*K(:,j);
    end %j    
    K(:,i) = F(x_it,y_it);
end %i
end

function [w,z] = getWZ(y,h,K,BW,BZ)
s = length(BW);
w = y; z = y;
for i = 1:s
   w = w + h*BW(i)*K(:,i);
   z = z + h*BZ(i)*K(:,i);
end
end

function dw = ex2(~,w) % generally expect 2 arguments for solvers (IV first)
% for ode45 etc... must return a column vector.
dw = nan(2,1);
dw(1) = w(2);
dw(2) = -0.25*(4*w(2) + 17*w(1));
end



