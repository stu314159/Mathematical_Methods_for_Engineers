%% Lecture 31 Example
clear
clc
close 'all'

%% Parameters
N = 50; % number of modes
c = 1; % radius of the circle
a_sq = 1.0; % "stiffness" parameter
a = sqrt(a_sq);

example = 2;
%example = [1 | 2]
switch example
    case 1        
        f = @(r) ex1(r,c); % initial position
        g = @(r) 0.*r; % initial velocity
    case 2
        b = 0.2;
        f = @(r) 0.*r;
        g = @(r) ex2(r,b);
    otherwise
        error('Unexpected example number!!');
end
%% Get eigenvalues
k = besselzero(0,N,1);
nu = k./c;

R = @(r,n) besselj(0,nu(n).*r);

U = @(r,t) 0;

for n = 1:N
   % compute An
   an = integral(@(r) r.*f(r).*R(r,n),0,c)./...
       integral(@(r) r.*(R(r,n).^2),0,c);
   % compute Bn
   bn = integral(@(r) r.*g(r).*R(r,n),0,c)./...
       (a.*nu(n).*integral(@(r) r.*R(r,n).^2,0,c));
   
   % add the term to our solution
   U = @(r,t) U(r,t) + (an*cos(a*nu(n).*t) + ...
       bn*sin(a*nu(n).*t)).*R(r,n);
end

%% Plot the solution
NR = 20;
NTHETA = 20;
Tmax = 10;
NT = 50;
R = linspace(0,c,NR);
THETA = linspace(0,2*pi,NTHETA);
[RR,TT] = meshgrid(R,THETA);
XX = RR.*cos(TT);
YY = RR.*sin(TT);

% T = linspace(0,Tmax,NT);
% 1.5
% for t = 1:NT
%    UUp = U(RR,T(t));
%    surf(XX,YY,UUp,'facecolor','none');
%    title_str = sprintf('Lecture 31 example, t = %g \n',T(t));
%    title(title_str,'fontsize',16,'fontweight','bold');
%    xlabel('X','fontsize',14,'fontweight','bold');
%    ylabel('Y','fontsize',14,1.0'fontweight','bold');
%    zlabel('U','fontsize',14,'fontweight','bold');
%    set(gca,'fontsize',12,'fontweight','bold');
%    axis([-(c+.1) c+.1 -(c+.1) c+.1 -2 2]);
%    pause(0.5*Tmax/(NT-1));  
%     
% end

%% fixed plots for inclusion in the book
t_i = 3.5;
surf(XX,YY,U(RR,t_i),'faceColor','none');
title_str = sprintf('Lecture 31 example, t = %g \n',t_i);
title(title_str,'fontsize',16,'fontweight','bold');
xlabel('X','fontsize',14,'fontweight','bold');
ylabel('Y','fontsize',14,'fontweight','bold');
zlabel('U','fontsize',14,'fontweight','bold');
set(gca,'fontsize',12,'fontweight','bold');
axis([-(c+.1) c+.1 -(c+.1) c+.1 -2 2]);


%% Local functions
function y = ex1(r,c)
[m,n] = size(r);
y = nan(m,n);
for i = 1:length(r)
   if (r(i) < c/3)
       y(i) = 1;
   else
       y(i) = 0;
   end    
end
end

function y = ex2(r,b)
[m,n] = size(r);
y = nan(m,n);
vo = 10;
for i = 1:length(r)
   if (r(i) < b)
       y(i) = -vo;
   else
       y(i) = 0;
   end    
end
end