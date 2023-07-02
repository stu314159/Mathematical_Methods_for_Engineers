%% Lecture 32 example 2
clear
clc
close 'all'

%% Case 2
R = 1; Z = 1;
g = @(z) 1-z; % temperature boundary condition
N = 150;

c = nan(N,1);
u = @(r,z) 0; % initialize the series
for n = 1:N
    c(n) = (1./besseli(0,n*pi)).*...
        integral(@(z) g(z).*sin(n*pi*z),0,Z)./...
        integral(@(z) sin(n*pi*z).*sin(n*pi*z),0,Z);
    
    % update the series with the next term
    u = @(r,z) u(r,z) + ...
        c(n)*besseli(0,n*pi*r).*sin(n*pi*z);
    
end

%% make surface plot
Rv = linspace(0,R,100);
Zv = linspace(0,Z,200);

[RR,ZZ] = meshgrid(Rv,Zv);
UU = u(RR,ZZ);

figure(1)
surf(Rv,Zv,UU,'edgecolor','none');
title('Laplacian in a Cylinder: Case 2',...
    'fontsize',18,'fontweight','bold');
xlabel('R','fontsize',16,'fontweight','bold');
ylabel('Z','fontsize',16,'fontweight','bold');
zlabel('U','fontsize',16,'fontweight','bold');
view([-63 15]);