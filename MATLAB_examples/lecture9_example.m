%% Lecture9_example_r2.m
clear
clc
close 'all'

% this version will try and do without symbolic series

%% Parameters
n = 25;
xMin = 0; xMax = 5;

%% Call a local Function to Generate Partial Series Solution
y = generate_series(n);

%% Plot the Power Series
figure(1)
fplot(y,[xMin, xMax],'linewidth',2);
title_str = sprintf('Lecture 9 Series n = %d',n);
title(title_str,'fontsize',18,'fontweight','bold');
xlabel('X','fontsize',16,'fontweight','bold');
ylabel('Y(X)','fontsize',16,'fontweight','bold');
set(gca,'fontsize',12,'fontweight','bold');
grid on

%% Numeric solution with ODE45
y0 = [5,1]; % y(0) and y'(0)
tSpan = [0,5]; % think of "x" as time ("t")
[t,y_num] = ode45(@lec9_example,tSpan,y0);

figure(2)
plot(t,y_num(:,1),'linewidth',2);
title('Numeric Solution','fontsize',18,'fontweight','bold');
xlabel('X','fontsize',16,'fontweight','bold');
ylabel('Y(X)','fontsize',16,'fontweight','bold');
set(gca,'fontsize',12,'fontweight','bold');
grid on

figure(3)
fplot(y,[xMin,xMax],'-b','linewidth',2);
hold on
plot(t,y_num(:,1),'--r','linewidth',2);
hold off
title_str = sprintf('Power Series Comparison n = %d',n);
title(title_str,'fontsize',18,'fontweight','bold');
xlabel('X','fontsize',16,'fontweight','bold');
ylabel('Y(X)','fontsize',16,'fontweight','bold');
set(gca,'fontsize',12,'fontweight','bold');
legend('Power Series','ODE45');
grid on

%% Compute difference between y_num and y
sol_rel_err = (y_num(:,1) - y(t))./y(t);
figure(4)
plot(t,sol_rel_err,'-c','LineWidth',3);
title('Relative Error','fontsize',14,'fontweight','bold');
xlabel('X','FontSize',12,'FontWeight','bold');
ylabel('Relative Error','FontSize',12,'FontWeight','bold');
sol_rel_err_norm = norm(sol_rel_err,2);
fprintf('The 2-norm of the relative error = %g',...
    sol_rel_err_norm);

%% Local function to generate the series
function  y = generate_series(n)

C1 = nan(1,n);
C2 = nan(1,n);

y1 = @(t) 0;
c1_0 = 5;
C1(1) = 0;
C1(2) = c1_0/2; % c2 = c0/2

k = 1;
C1(k+2) = (C1(k)+c1_0)/((k+1)*(k+2));

% compute the coefficients
for k = 2:(n-3)
    C1(k+2) = (C1(k) + C1(k-1))/((k+1)*(k+2));
end

% construct y1 from the coefficients
y1 = @(t) y1(5) + c1_0;
for k = 1:(n-1)
    y1 = @(t) y1(t) + C1(k)*t.^k;
end
% for y2(t), c0 = 0, c1 = 1

% set coefficients for y2(t)
c2_0 = 0;
C2(1) = 1;
C2(2) = c2_0/2; % again, c2 = c0/2

k = 1;
C2(k+2) = (C2(k) + c2_0)/((k+1)*(k+2));
% compute the coefficients
for k = 2:(n-3)
   C2(k+2)= (C2(k) + C2(k-1))/((k+1)*(k+2)); 
end

% construct y2 from the coefficients
y2 = @(t) c2_0;
for k = 1:(n-1)
    y2 = @(t) y2(t) + C2(k)*t.^k;
end

y = @(t) y1(t) + y2(t);


end

%%  Local function for generation of numeric solution with ODE45
function dw_dt = lec9_example(t,w)

dw_dt = zeros(2,1);
dw_dt(1) = w(2);
dw_dt(2) = (1+t).*w(1);

end
