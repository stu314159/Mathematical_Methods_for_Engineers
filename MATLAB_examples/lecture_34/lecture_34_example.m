%% Lecture 34 - Laplacian in Spherical Coordinates example
% Steady state temperature profile for a sphere.
clear
clc
close 'all'
%% Set Parameters
R = 2;
N = 4;
f = @(theta) ex1(theta);

%% Construct the Solution
c = nan(N,1);
u = @(r,theta) 0; % initialize the series

% start for n = 0 
n = 0;
c0 = integral(@(th) f(th).*...
    legendreP(n,cos(th)).*sin(th),0,pi)./...
    integral(@(th) (R^n)*...
    (legendreP(n,cos(th)).^2).*sin(th),0,pi);

u = @(r,theta) u(r,theta) + c0*(r.^0).*legendreP(n,cos(theta));    

for n = 1:N 
    % get the next coeficient
    c(n) = integral(@(th) f(th).*...
        legendreP(n,cos(th)).*sin(th),0,pi)./...
        integral(@(th) (R^n)*...
        (legendreP(n,cos(th)).^2).*sin(th),0,pi);
    
    % update the approximation
    u = @(r,theta) u(r,theta) + ...
        c(n)*(r.^n).*legendreP(n,cos(theta)); 
end

%% Process Result for Plotting 
Nx = 50;
Xv = linspace(-R,R,Nx);
Yv = linspace(-R,R,Nx);
Zv = linspace(-R,R,Nx);
dx = Xv(2)-Xv(1); % need this for VTK file
[XX,YY,ZZ] = meshgrid(Xv,Yv,Zv);

RR = sqrt(XX.^2+YY.^2+ZZ.^2);
PP = acos(ZZ./RR);
UU = u(RR,PP); 
% set region outside the sphere to nan
UU(RR>R) = nan;

%% Write the data to a file
filename = 'solution.vtk';
dataname = 'U';
origin = [-R -R -R]; 
spacing = [dx dx dx];
save_scalarStructuredPoints3D_VTK_binary(filename,...
    dataname,UU,origin,spacing);

%% Local functions
function u = ex1(theta)
[n,m] = size(theta);
u = nan(n,m);
ind_a = theta<=(pi/2);
ind_b = theta>(pi/2);
u(ind_a) = 10*((pi/2)-theta(ind_a));
u(ind_b) = 0;
end

function save_scalarStructuredPoints3D_VTK_binary(filename,...
    dataname,data_set,origin,spacing)

[nx,ny,nz]=size(data_set);

% open the file
fid = fopen(filename,'w');

% ASCII file header
fprintf(fid,'# vtk DataFile Version 3.0\n');
fprintf(fid,'VTK from Matlab\n');
fprintf(fid,'BINARY\n\n');
fprintf(fid,'DATASET STRUCTURED_POINTS\n');
fprintf(fid,'DIMENSIONS %d %d %d \n',nx,ny,nz);
fprintf(fid,'ORIGIN  %4.3f   %4.3f  %4.3f \n',...
    origin(1),origin(2),origin(3));
fprintf(fid,'SPACING %4.3f   %4.3f  %4.3f \n',...
    spacing(1),spacing(2),spacing(3));
fprintf(fid,'\n');
fprintf(fid,'POINT_DATA %d \n',nx*ny*nz);
fprintf(fid,strcat('SCALARS','\t ',dataname,' float ', '\n'));
fprintf(fid,'LOOKUP_TABLE default \n');
% write the data
fwrite(fid,reshape(data_set,1,nx*ny*nz),'float','b');
% close the file
fclose(fid);

end
