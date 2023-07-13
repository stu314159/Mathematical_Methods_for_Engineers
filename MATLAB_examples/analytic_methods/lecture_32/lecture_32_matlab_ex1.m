%% Lecture 32 example 1
clear
clc
close 'all'

%% Parameters
R = 2; Z = 4;
Uo = 5; % given temperature on top surface of the cylinder
N = 100;
k = besselzero(0,N,1); % get first n zeros of Jo
alpha = k/R;

%% Construct Solution
A = nan(N,1);
u = @(r,z) 0; % initialize the series
for n = 1:N
    A(n) = (Uo/(sinh(Z*alpha(n)))).*...
        integral(@(r) besselj(0,alpha(n)*r).*r,0,R)./...
        integral(@(r) besselj(0,alpha(n)*r).*...
        besselj(0,alpha(n)*r).*r,0,R);
    
    % update the series with the next term
    u = @(r,z) u(r,z) + ...
        A(n)*besselj(0,alpha(n)*r).*sinh(z*alpha(n));    
end


%% Plot results
Rv = linspace(0,R,100);
Zv = linspace(0,Z,200);

[RR,ZZ] = meshgrid(Rv,Zv);
UU = u(RR,ZZ);

figure(1)
surf(Rv,Zv,UU,'edgecolor','none');
title('Laplacian in a Cylinder: Case 1',...
    'fontsize',18,'fontweight','bold');
xlabel('R','fontsize',16,'fontweight','bold');
ylabel('Z','fontsize',16,'fontweight','bold');
zlabel('U','fontsize',16,'fontweight','bold');
view([65 10]);