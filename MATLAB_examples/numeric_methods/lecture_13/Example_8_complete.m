%% Example 8 - linear least squares
clear
clc
close 'all'

%% Input Data
I = [300 300 350 400 400 500 500 650 650]';
d = [22 26 27 30 34 33 33.5 37 42]';

%% Plot the data
figure(1)
plot(I,d,'sk','markersize',10,...
    'linewidth',2);
title('Raw Data','fontsize',18,...
    'fontweight','bold');
xlabel('Current Flow (nA)',...
    'fontsize',16,'fontweight','bold');
ylabel('Fiber Diameter (\mum)',...
    'fontsize',16,'fontweight','bold');
grid on
set(gca,'fontsize',12,'fontweight','bold');
axis([250 675 20 45]);

%% Linear Least Squares
fprintf('\n\n Linear Least Squares \n');
X = [ I.^0 I.^1];
b = d;

c = (X'*X)\(X'*b);

fprintf('b = %g \n',c(1));
fprintf('m = %g \n',c(2));

linEst = @(x) c(1) + c(2)*x;
% plot the fit line
figure(2)
plot(I,d,'sk',...
    I,linEst(I),'-r','linewidth',3,...
    'markersize',10);
title('Linear Fit','fontsize',18,...
    'fontweight','bold');
xlabel('Current Flow (nA)','fontsize',16,...
    'fontweight','bold');
ylabel('Fiber Diameter (\mum)','fontsize',16,...
    'fontweight','bold');
grid on
set(gca,'fontsize',12,'fontweight','bold');
axis([250 675 20 45]);


%% Least squares - second order

fprintf('\n\n Quadradic Least Squares \n');
X = [I.^0 I.^1 I.^2];
b = d;

c = (X'*X)\(X'*b);

fprintf('a = %g \n',c(1));
fprintf('b = %g \n',c(2));
fprintf('c = %g \n',c(3));

quadEst = @(x) c(1) + c(2)*x + c(3)*x.^2;
% plot the fit line
figure(3)
plot(I,d,'sk',...
    I,linEst(I),'-r',...
    I,quadEst(I),'-g','linewidth',3,'markersize',10);
title('Second-Order Fit','fontsize',18,'fontweight','bold');
xlabel('Current Flow (nA)','fontsize',16,'fontweight','bold');
ylabel('Fiber Diameter (\mum)','fontsize',16,'fontweight','bold');
grid on
set(gca,'fontsize',12,'fontweight','bold');
axis([250 675 20 45]);

%% Least Squares - (problem 6.9)

fprintf('\n\n Least Squares - Model Fit \n');
X = [I.^0 I.^0.5];
b = d;
c = (X'*X)\(X'*b);
fprintf('a = %g \n',c(1));
fprintf('b = %g \n',c(2));

est3 = @(x) c(1) + c(2)*x.^0.5;

% plot the fit line
figure(4)
plot(I,d,'sk',...
    I,linEst(I),'-r',...
    I,quadEst(I),'-g',...
    I,est3(I),'-c','linewidth',3,'markersize',10);
title('Model Fit','fontsize',18,'fontweight','bold');
xlabel('Current Flow (nA)','fontsize',16,'fontweight','bold');
ylabel('Fiber Diameter (\mum)','fontsize',16,'fontweight','bold');
grid on
set(gca,'fontsize',12,'fontweight','bold');
legend('Data','Linear Fit','Quad Fit','Model Fit');
axis([250 675 20 45]);