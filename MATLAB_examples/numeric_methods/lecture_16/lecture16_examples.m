%% Lecture 16 Examples
clear
clc
close 'all'

%% Example 1
base_year = 1981;
year = [1981,1984, 1989, 1993, 1997,...
    2000, 2001, 2003, 2004, 2010];


pct_w_comp = [0.5, 8.2, 15, 22.9, 36.6, 51,...
    56.3, 61.8, 65, 76.7];

figure(1)
plot(year,pct_w_comp,'ro','markersize',10);
grid on
title('Percent Household with Comuters vs. Year',...
    'fontsize',14,'fontweight','bold');
xlabel('Year','fontsize',12,'fontweight','bold');
ylabel('% Households with Computers','fontsize',12,...
    'fontweight','bold');

year = year-base_year;
m = 3;
p1 = polyfit(year,pct_w_comp,m);

x1 = 2008 - base_year; x2 = 2013 - base_year;

fprintf('Estimated percent ownership in %d is %g percent.\n',...
    x1+base_year, polyval(p1,x1));
fprintf('Estimated percent ownership in %d is %g percent.\n',...
    x2+base_year, polyval(p1,x2));

yearMax = max(year); yearMin = min(year);
N = 100;
Years = linspace(yearMin,yearMax,N);
figure(2)
plot(year+base_year,pct_w_comp,'ro',...
    Years+base_year,polyval(p1,Years),'-b',...
    'markersize',10);
grid on
title('Percent Household with Comuters vs. Year',...
    'fontsize',14,'fontweight','bold');
xlabel('Year','fontsize',12,'fontweight','bold');
ylabel('% Households with Computers','fontsize',12,...
    'fontweight','bold');

%% Example 2
x = [8, 11, 15, 18, 22];
y = [5, 9, 10, 8, 7];

method = 'pchip';

f_int = @(xi) interp1(x,y,xi,method);

xMin = min(x); xMax = max(x); Nx = 1000;
X_space = linspace(xMin, xMax, Nx);

figure(3)
plot(x,y,'ro',...
    X_space,f_int(X_space),'-b',...
    'linewidth',3)
title('Interpolation Methods','fontsize',16,...
    'fontweight','bold');
xlabel('X','fontsize',14,'fontweight','bold');
ylabel('Y','fontsize',14,'fontweight','bold');
set(gca,'fontsize',12,'fontweight','bold');
axis([xMin xMax 0.5*min(y) 1.25*max(y)]);

%% Example 3
x = [1200, 1500, 2000, 2500, 3000,...
    3250, 3500, 3750,4000, 4400];
y = [65, 130, 185, 225, 255, 266,...
    275, 272, 260, 230];

method = 'pchip';

f_int = @(xi) interp1(x,y,xi,method);

xMin = min(x); xMax = max(x); Nx = 1000;
X_space = linspace(xMin, xMax, Nx);

figure(4)
plot(x,y,'ro',...
    X_space,f_int(X_space),'-b',...
    'linewidth',3)
title('MATLAB Interpolation','fontsize',16,...
    'fontweight','bold');
xlabel('Engine Speed [RPM]','fontsize',14,'fontweight','bold');
ylabel('Power [hp]','fontsize',14,'fontweight','bold');
grid on;
set(gca,'fontsize',12,'fontweight','bold');
axis([xMin xMax 0.5*min(y) 1.25*max(y)]);

x1 = 2300;
fprintf('Estimated hp at %d rpm: %g \n',...
    x1, f_int(x1));
x2 = 3650;
fprintf('Estimated hp at %d rpm: %g \n',...
    x2,f_int(x2));

%% Example 4
% load data from FEM analysis
load('fem_data.mat');

figure(5)
scatter3(gcoord(:,1),gcoord(:,2),T,15,T,'filled');

TempField = scatteredInterpolant(gcoord(:,1),gcoord(:,2),T);

xMin = min(gcoord(:,1)); xMax = max(gcoord(:,1));
yMin = min(gcoord(:,2)); yMax = max(gcoord(:,2));

Nx = 200;
Ny = 200;
X = linspace(xMin,xMax,Nx);
Y = linspace(yMin,yMax,Ny);
[XX,YY] = meshgrid(X,Y);
figure(6)
surf(XX,YY,TempField(XX,YY),'edgecolor','none');
title('Temperature Field','fontsize',16,...
    'fontweight','bold')
xlabel('X','fontsize',14,'fontweight','bold');
ylabel('Y','FontSize',14, 'FontWeight','bold');
zlabel('T ^{\circ}C','FontSize',14,'FontWeight','bold');
axis([xMin xMax yMin yMax 400 700]);
