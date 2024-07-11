%% Lecture_30_matlab.m

clear
clc
close 'all'

%% Parameters
c = 2;
N = 50;

fx_pick = 1;
%[1 | 2]
switch fx_pick
    case 1
        f = @(x) x;
    case 2
        f = @(x) ex1(x);
    otherwise
        error('Invalid case!');    
end

%% solve the problem

Ao = (1./(2*pi))*integral(@(theta) f(theta),0,2*pi);

U = @(r,theta) Ao;
reltol = 1e-15;
for n = 1:N
   
    An = (1/((c.^n)*pi))*integral(@(theta) f(theta).*cos(n.*theta),...
        0,2*pi,'RelTol',reltol); %<-- set alternate relative tolerance 
    Bn = (1/((c.^n)*pi))*integral(@(theta) f(theta).*sin(n.*theta),...
        0,2*pi,'RelTol',reltol);  %<--  ditto.
    
   
    U = @(r,theta) U(r,theta)+ (r.^n).*(An*cos(n*theta)+Bn*sin(n*theta));
    
end

%% Make a Plot
NR = 100;
NT = 100;
R = linspace(0,c,NR);
THETA = linspace(0,2*pi,NT);
[RR,TT] = meshgrid(R,THETA);
UUp = U(RR,TT);

% plot in cartesian coordinates
XX = RR.*cos(TT); % get cartesian coordinate equivalents
YY = RR.*sin(TT);

figure(1)
surf(XX,YY,UUp,'edgecolor','none');
colormap('jet');%<-- consider alternate colormaps
c = colorbar;%<-- add a colorbar
c.Label.String = 'Temperature'; %<-- give colorbar a label
title('Lecture 30 Example','fontsize',16,'fontweight','bold');
xlabel('X','fontsize',14,'fontweight','bold');
ylabel('Y','fontsize',14,'fontweight','bold');
zlabel('U','fontsize',14,'fontweight','bold');
grid on
set(gca,'fontsize',12,'fontweight','bold');

%% Local functions
function y = ex1(theta)
[m,n] = size(theta);
y = nan(m,n);
for i = 1:length(theta)
    if(theta(i)>= 0) && (theta(i)< pi/2)
        y(i) = 1;
    else
        y(i) = 0;
    end
end    
end

