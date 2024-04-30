%% Assignment 5 Solution script
clear
clc
close 'all'

%% Problem 1
fprintf('\n\n Problem #1 \n\n');
% part a)
fprintf('Part a) \n');

x = [1 3 4 6 9 12 14]';% degrees C, given
y = [2 4 5 6 7 9 11]';% atmospheres
fprintf('\n\n Problem #1 solution using LinReg\n\n');
[a,Er] = LinReg(x,y);
fprintf('Coefficients: \n');
disp(a);
snapnow

fprintf('Square Residual (error): %g \n',Er);
snapnow

%% Problem 2

fprintf('\n\n Problem #2 \n\n');

T = [595 623 761 849 989 1076 1146 1202 1382 1445 1562]';
k = (1e-20)*[2.12 3.12 14.4 30.6 80.3 131 186 240 489 604 868]';

F1 = @(x) x.^0;
F2 = @(x) log(x);
F3 = @(x) 1./x;

% get the fit.
C = NonLinCombFit(F1,F2,F3,T,log(k));

% undo the transformation
nonLinEst = @(x) exp(C(1) + C(2)*log(x) + C(3)./x);

fprintf('Coefficients: \n');
disp(C);
snapnow

figure(2)
plot(T,k,'sr',...
    T,nonLinEst(T),'-c','linewidth',3,'markersize',12);
title('Problem #2 Solution','fontsize',18,'fontweight','bold');
xlabel('T (K)','fontsize',16,'fontweight','bold');
ylabel('k (m^3/s) ','fontsize',16,'fontweight','bold');
grid on
set(gca,'fontsize',12,'fontweight','bold');

% A = exp(C(1)); -E/R = C(3)--> E = -R*C(3)
fprintf('A = %g \n',exp(C(1)));
R = 8.314; % J/mole/K
E = -R*C(3);
fprintf('E = %g J/mole\n',E);

%% Problem #3 (check)
x = [0.2 0.5 1 2 3]';
y = [3 2 1.4 1 0.6]';

 % transform 1/y = mx + b
 p = 1./y;
 X = [x.^0 x.^1];
 
 c = (X'*X)\(X'*p);
 
 b = c(1);
 m = c(2);
 
 est = @(x) 1./(m*x + b);
 
 figure(3)
 plot(x,y,'sr',...
     x, est(x),'-b','linewidth',3,'markersize',8)
 grid on
 title('Problem #3 Check','fontsize',16,'fontweight','bold');
     




%% Problem #4
fprintf('\n\n Problem #4 \n\n');

T = [20 100 180 260 340 420 500]'; % degrees C
R = [500 676 870 1060 1205 1410 1565]';% ohms resistance

To = 20; % degrees C


% Part a)
fprintf('\n Part a) \n');
% R = Ro*(1 + alpha*(T - To)] <-- find Ro and alpha
% R = Ro + Ro*alpha*(T-To)
% send to LinReg to get Ro and alpha
% C(1) will be Ro, C(2) will be Ro*alpha

[C,Er] = LinReg((T-To),R);

fprintf('Coefficients: \n');
disp(C);
fprintf('Error estimate: %g \n',Er);

Ro = C(1);
alpha = C(2)/Ro;

fprintf('Ro = %g \n',Ro);
fprintf('alpha = %g \n',alpha);

linEst = @(T) Ro*(1 + alpha*(T - To)); %<- note T in this equation is different than T outside.
figure(3)
plot(T,R,'sr',...
    T,linEst(T),'-c','linewidth',3,'markersize',12);
title('Problem #4 Solution','fontsize',18,'fontweight','bold');
xlabel('T (C)','fontsize',16,'fontweight','bold');
ylabel('Resistance (ohms) ','fontsize',16,'fontweight','bold');
grid on
set(gca,'fontsize',12,'fontweight','bold');

%% Problem #5

fprintf('\n\n Problem #5 \n\n');

F = 1e3*[24.6 29.3 31.5 33.3 34.8 35.7 36.6 37.5 38.8 39.6 40.4]'; % N
L = 1e-3*[12.58 12.82 12.91 12.95 13.05 13.21 13.35 13.49 ...
    14.08 14.21 14.48]'; % m

Ao = 1.25e-4; %initial cross sectional area, m^2
Lo = 0.0125; % inital length, m

stress = F.*L./(Ao*Lo);
strain = log(L/Lo);

X = [ log(strain).^0 log(strain)];
b = log(stress);

C = (X'*X)\(X'*b);

% C(1) = log(k)
% C(2) = m

k = exp(C(1));
m = C(2);

est = @(x) k*(x.^m);

fprintf('K = %g \n',k);
fprintf('m = %g \n',m);

figure(5)
plot(strain,stress,'sr',...
    strain,est(strain),'-c','linewidth',3,'markersize',12);
title('Problem #5 Solution ','fontsize',18,'fontweight','bold');
xlabel('Strain','fontsize',16,'fontweight','bold');
ylabel('Stress (Pa) ','fontsize',16,'fontweight','bold');
grid on
set(gca,'fontsize',12,'fontweight','bold');


%% Local Functions
function [a,Er] = LinReg(x,y)
% function [a,Er] = linReg(x,y) performs linear regression to fit data x
% and y.  
%
% Inputs:
% x - vector with independent values
% y - vector with dependent values
% 
% Outputs:
% a - vector of coefficients to linear fit to x,y
% Er - scalar representing the squared residual

% I absolutely refuse to use the implementation provided in the text.
[N,M] = size(x);
N_len = length(x);
if N ~= N_len && M == N_len % I want x to be a Nx1 vector
    x = x';
end

[N,M] = size(y);
N_len = length(y);
if N ~= N_len && M == N_len % I want y to be a Nx1 vector
    y = y';
end

X = [x.^0 x.^1];
a = (X'*X)\(X'*y);

r = y - (a(2)*x + a(1));
Er = r'*r;
end

function C = NonLinCombFit(F1,F2,F3,x,y)

% I absolutely refuse to use the implementation provided in the text.
[N,M] = size(x);
N_len = length(x);
if N ~= N_len && M == N_len % I want x to be a Nx1 vector
    x = x';
end

[N,M] = size(y);
N_len = length(y);
if N ~= N_len && M == N_len % I want y to be a Nx1 vector
    y = y';
end

X = [F1(x) F2(x) F3(x)];

C = (X'*X)\(X'*y);


end

