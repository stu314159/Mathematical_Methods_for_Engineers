%% Lecture 15 Demos

clear
clc
close 'all'

%% Data
% for example shown on pg 212
x = [1 4 7 10 13]';
y = [2 6 4 8 10]';

%% Plot the data
figure(1)
plot(x,y,'rs','linewidth',3,'markersize',10);
title('Sample Data','fontsize',18,'fontweight','bold');
xlabel('X','fontsize',16,'fontweight','bold');
ylabel('Y','fontsize',16,'fontweight','bold');
grid on
set(gca,'fontsize',12,'fontweight','bold');
axis([0 15 0 15]);

%% Polynomial Interpolation

N = length(x);
X = x.^(0:N);
c = (X'*X)\(X'*y);
nthInterp = @(x) (x.^(0:N))*c;

x_plt = linspace(min(x),max(x),1000);
x_plt = x_plt';

figure(2)
plot(x,y,'sr',...
    x_plt,nthInterp(x_plt),'-k','markersize',10,'linewidth',3);
title('Least Squares Interpolation','fontsize',18,'fontweight','bold');
xlabel('X','fontsize',16,'fontweight','bold');
ylabel('Y','fontsize',16,'fontweight','bold');
grid on
legend('Sample Data','Interpolant')
set(gca,'fontsize',12,'fontweight','bold');
axis([0 15 0 15]);


%% Data
% from example 6.4 on page 215
x = [1 2 4 5 7]';
y = [52 5 -5 -40 10]';

%% Plot the data
figure(3)
plot(x,y,'rs','linewidth',3,'markersize',10);
title('Sample Data','fontsize',18,'fontweight','bold');
xlabel('X','fontsize',16,'fontweight','bold');
ylabel('Y','fontsize',16,'fontweight','bold');
grid on
set(gca,'fontsize',12,'fontweight','bold');
axis([0 8 -45 55]);

%% Interpolation as the logical limit of data fitting
%%
% 
%  Having just done a section on data fitting, we can consider what happens
%  as we generate successively more detailed fits of the given data
% 

% linear
X = x.^(0:1);%<-- does this make sense?
c = (X'*X)\(X'*y);

linInterp = @(x) (x.^(0:1))*c; %<-- likewise, this?
figure(4)
plot(x,y,'sr',...
    x,linInterp(x),'-r','markersize',10,'linewidth',3);
title('Sample Data','fontsize',18,'fontweight','bold');
xlabel('X','fontsize',16,'fontweight','bold');
ylabel('Y','fontsize',16,'fontweight','bold');
grid on
set(gca,'fontsize',12,'fontweight','bold');
axis([0 8 -45 55]);

%% N-th order
%%
% 
%  We can make a successively more "accurate" approximation of the function
%  by using a higher-order linear data fit.  With 5 data points, we can
%  match it exactly with a 4th-order polynomial (e.g. 5 parameters)
% 

N = 5;
X = x.^(0:N);
c = (X'*X)\(X'*y);
nthInterp = @(x) (x.^(0:N))*c;
figure(5)
plot(x,y,'sr',...
    x,nthInterp(x),'-r','markersize',10,'linewidth',3);
title('Sample Data','fontsize',18,'fontweight','bold');
xlabel('X','fontsize',16,'fontweight','bold');
ylabel('Y','fontsize',16,'fontweight','bold');
grid on
set(gca,'fontsize',12,'fontweight','bold');
axis([0 8 -45 55]);

%%
% 
%  Thus if our goal was interpolation, we could use the 4th order data fit
%  as an interpolant.  The coefficients found (c) are the coefficients for
%  a 4th order polynomial that interpolates the data.
% 

%%
% 
%  This is problematic for the following reasons:
% 
%%
% 
% # Determining the coefficients is a little bit of trouble (solve system
% of equations
% # The system of equations that you solve gets less "well conditioned" as
% N increases.
% 

%% Interpolation Demo (part 2)
% Based on Example 6-3 on page 208 of the textbook

clear
clc
close 'all'

%% Load data
Strain = [0 0.4 0.8 1.2 1.6 2.0 2.4 2.8 3.2...
    3.6 4.0 4.4 4.8 5.2 5.6 6.0]';% dimensionless
Stress = [0 3.0 4.5 5.8 5.9 5.8 6.2 7.4 9.6...
    15.6 20.7 26.7 31.1 35.6 39.3 41.5]'; % MPa

strplt = linspace(min(Strain),max(Strain),1000);
strplt = strplt';

%% Make a Basic Plot of the Data
figure(1)
plot(Strain, Stress,'sr','markersize',10,'linewidth',3);
xlabel('Strain','fontsize',16,'fontweight','bold');
ylabel('Stress (MPa)','fontsize',16,'fontweight','bold');
title('Stress-Strain Curve','fontsize',18,'fontweight','bold');
grid on
set(gca,'fontsize',12,'fontweight','bold');

%% Linear Interpolation

X = Strain.^(0:1); %<-- less typing this way

c = (X'*X)\(X'*Stress);

linearInterp = @(x)  (x.^(0:1))*c; %<-- less typing and more general

figure(2)
plot(Strain, Stress,'sr',...
    Strain,linearInterp(Strain),'-r','markersize',10,'linewidth',3);
xlabel('Strain','fontsize',16,'fontweight','bold');
ylabel('Stress (MPa)','fontsize',16,'fontweight','bold');
title('Stress-Strain Curve','fontsize',18,'fontweight','bold');
grid on
set(gca,'fontsize',12,'fontweight','bold');

%% N-th order Interpolation (monomials)

N = 15;
X = Strain.^(0:N);
c = (X'*X)\(X'*Stress);
nthInterp = @(x) (x.^(0:N))*c;

figure(3)
plot(Strain, Stress,'sr',...
    strplt,nthInterp(strplt),'-r','markersize',10,'linewidth',3);
xlabel('Strain','fontsize',16,'fontweight','bold');
ylabel('Stress (MPa)','fontsize',16,'fontweight','bold');
title('Stress-Strain Curve','fontsize',18,'fontweight','bold');
grid on
set(gca,'fontsize',12,'fontweight','bold');

r = Stress - nthInterp(Strain);
fprintf('Residual Squared = %g \n',r'*r);

%% Lagrange Interpolation
F = genLagrangePolyInterp(Strain,Stress);



figure(4)
plot(Strain,Stress,'sr',...
    strplt,F(strplt),'-r','markersize',10,'linewidth',3);

xlabel('Strain','fontsize',16,'fontweight','bold');
ylabel('Stress (MPa)','fontsize',16,'fontweight','bold');
title('Stress-Strain Curve','fontsize',18,'fontweight','bold');
grid on
set(gca,'fontsize',12,'fontweight','bold');
r = Stress - F(Strain);
fprintf('Residual Squared = %g \n',r'*r);

%% Local function for Lagrange Polynomial
function F = genLagrangePolyInterp(X,Y)
% function F = genLagrangePoly(X,Y) generates a Lagrange Polynomial that
% may be used to interpolate a function
% Inputs
% X = x-values of a function
% Y = f(X) for some function
%
% Outputs
% F - a function handle with the Lagrange interpolant

n = length(X); 

F = @(x) 0; %<-- start with zero

for i = 1:n
    L = @(x) Y(i); %<-- initialize the Lagrange Function
    for j = 1:n
        if j ~= i
            L = @(x) L(x).*((x - X(j))./(X(i) - X(j)));
        end
    end
    F = @(x) F(x)+L(x);
end
end