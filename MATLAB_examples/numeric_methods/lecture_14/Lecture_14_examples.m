%% Example 9 - more least squares regression

clear
clc
close 'all'

%% Problem 1
fprintf('\nProblem 6.6 \n');

% data
x = [1 2 3 5 8]';
y = [0.5 1.9 2.2 3 3.5]';

% plot the data
figure(1)
plot(x,y,'sk','markersize',10,'linewidth',2);
title('Least Squares Problem 6.6 ','fontsize',18,'fontweight','bold');
xlabel('X','fontsize',16,'fontweight','bold');
ylabel('Y','fontsize',16,'fontweight','bold');
grid on
set(gca,'fontsize',12,'fontweight','bold');
axis([0 9 0 4]);

% linearized estimator
%X = 
%b = 

C = (X'*X)\(X'*b);
fprintf('m = %g \n',C(2));
fprintf('b = %g \n',C(1));

%remember to "undo" the linearizing transformation
%est1 = @(x) 

figure(1)
plot(x,y,'sk',...
    x,est1(x),'-r','markersize',10,'linewidth',2);
title('Least Squares Problem 6.6 ','fontsize',18,'fontweight','bold');
xlabel('X','fontsize',16,'fontweight','bold');
ylabel('Y','fontsize',16,'fontweight','bold');
grid on
set(gca,'fontsize',12,'fontweight','bold');
axis([0 9 0 4]);

%% Problem #2
fprintf('\n\nProblem 6.5 \n');

% data
x = [-2 -1 0 1 2]';
y = [1.5 3.2 4.5 3.4 2]';

% plot the data
figure(2)
plot(x,y,'sk','markersize',10,'linewidth',2);
title('Least Squares Problem 6.6 ','fontsize',18,'fontweight','bold');
xlabel('X','fontsize',16,'fontweight','bold');
ylabel('Y','fontsize',16,'fontweight','bold');
grid on
set(gca,'fontsize',12,'fontweight','bold');
axis([-2.5 2.5 1.5 5]);

% linearized estimator
%X = 
%b = 

C = (X'*X)\(X'*b);

% remember to "undo" the linearizing transformation
%a = 
%b = 

fprintf('a = %g \n',a);
fprintf('b = %g \n',b);

%est2 = @(x) 
% plot again
figure(2)
plot(x,y,'sk',...
    x,est2(x),'-r','markersize',10,'linewidth',2);
title('Least Squares Problem 6.6 ','fontsize',18,'fontweight','bold');
xlabel('X','fontsize',16,'fontweight','bold');
ylabel('Y','fontsize',16,'fontweight','bold');
grid on
set(gca,'fontsize',12,'fontweight','bold');
axis([-2.5 2.5 1.5 5]);

%% Problem #3
fprintf('\n\nProblem 6.8 \n');

% data
T = [-40 -20 0 20 40]';
W = [0.0012 0.002 0.0032 0.006 0.0118]';

% plot the data
figure(3)
plot(T,W,'sk','markersize',10,'linewidth',2);
title('Least Squares Problem 6.8 ','fontsize',18,'fontweight','bold');
xlabel('Temperature ^{\circ}C','fontsize',16,'fontweight','bold');
ylabel('Water Solubility (% wt.)','fontsize',16,'fontweight','bold');
grid on
set(gca,'fontsize',12,'fontweight','bold');
axis([-45 45 0 0.0135]);

%X = 
%b = ;

C = (X'*X)\(X'*b);


%b = ;%<-- note: over-writing other value of b
%m = 
fprintf('b = %g \n',b);
fprintf('m = %g \n',m);

%est3 = @(x) 

% plot the data again
figure(3)
plot(T,W,'sk',...
    T,est3(T),'-r','markersize',10,'linewidth',2);
title('Least Squares Problem 6.8 ','fontsize',18,'fontweight','bold');
xlabel('Temperature ^{\circ}C','fontsize',16,'fontweight','bold');
ylabel('Water Solubility (% wt.)','fontsize',16,'fontweight','bold');
grid on
set(gca,'fontsize',12,'fontweight','bold');
axis([-45 45 0 0.0135]);

fprintf('Solubility at 10 degrees C = %g \n',est3(10));

