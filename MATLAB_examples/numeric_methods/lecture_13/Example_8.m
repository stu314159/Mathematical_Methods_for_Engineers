%% Example 8 - linear least squares
clear
clc
close 'all'

%% Input Data
I = [300 300 350 400 400 500 500 650 650]';
d = [22 26 27 30 34 33 33.5 37 42]';

%% Plot the data
figure(1)
plot(I,d,'sk','markersize',10,'linewidth',2);
title('Raw Data','fontsize',18,'fontweight','bold');
xlabel('Current Flow (nA)','fontsize',16,'fontweight','bold');
ylabel('Fiber Diameter (\mum)','fontsize',16,'fontweight','bold');
grid on
set(gca,'fontsize',12,'fontweight','bold');
axis([250 675 20 45]);

%% Linear Least Squares
fprintf('\n\n Linear Least Squares \n');
%X = 
%b = 

%c = 

fprintf('b = %g \n',c(1));
fprintf('m = %g \n',c(2));

%linEst =
% plot the fit line
figure(2)
plot(I,d,'sk',...
    I,linEst(I),'-r','linewidth',3,'markersize',10);
title('Raw Data','fontsize',18,'fontweight','bold');
xlabel('Current Flow (nA)','fontsize',16,'fontweight','bold');
ylabel('Fiber Diameter (\mum)','fontsize',16,'fontweight','bold');
grid on
set(gca,'fontsize',12,'fontweight','bold');
axis([250 675 20 45]);


%% Least squares - second order

fprintf('\n\n Quadradic Least Squares \n');
%X = 
%b = 

%c = 

fprintf('a = %g \n',c(1));
fprintf('b = %g \n',c(2));
fprintf('c = %g \n',c(3));

%quadEst = 
% plot the fit line
figure(3)
plot(I,d,'sk',...
    I,linEst(I),'-r',...
    I,quadEst(I),'-g','linewidth',3,'markersize',10);
title('Raw Data','fontsize',18,'fontweight','bold');
xlabel('Current Flow (nA)','fontsize',16,'fontweight','bold');
ylabel('Fiber Diameter (\mum)','fontsize',16,'fontweight','bold');
grid on
set(gca,'fontsize',12,'fontweight','bold');
axis([250 675 20 45]);

%% Least Squares - (problem 6.9)

fprintf('\n\n Least Squares - non-monomial \n');
%X = 
%b = 
%c = 
fprintf('a = %g \n',c(1));
fprintf('b = %g \n',c(2));

%est3 = 

% plot the fit line
figure(4)
plot(I,d,'sk',...
    I,linEst(I),'-r',...
    I,quadEst(I),'-g',...
    I,est3(I),'-c','linewidth',3,'markersize',10);
title('Raw Data','fontsize',18,'fontweight','bold');
xlabel('Current Flow (nA)','fontsize',16,'fontweight','bold');
ylabel('Fiber Diameter (\mum)','fontsize',16,'fontweight','bold');
grid on
set(gca,'fontsize',12,'fontweight','bold');
axis([250 675 20 45]);