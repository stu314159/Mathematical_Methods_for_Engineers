%% Lecture 25 MATLAB - Heat Eqn
clear
clc
close 'all'

%% Heat Equation BVP & Solution
L = 10; % cm
alpha_sq = 1.752; % cm^2/s, thermal diffusivity of silver.

N = 50;

F = @(x,n) sin(n.*pi.*x./L);
G = @(t,n) exp(-((n.*pi./L).^2)*alpha_sq.*t);

f = @(x) 4 - 0.8*abs(x - 5);

u = @(x,t) 0;
for n = 1:N
    % essentially doing the sine-series half-wave expansion
    % compute the coefficient
    cn = (2/L)*integral(@(x) f(x).*F(x,n),0,L);
    
    % add the term to the series solution
    u = @(x,t) u(x,t) + cn*F(x,n).*G(t,n);
end


%% Plot the result for fixed times
% make a discrete X-axis
Nx = 1000;
X = linspace(0,L,Nx);
figure(1)
plot(X,u(X,0),'-ob',...
    X,u(X,1),'-.g',...
    X,u(X,10),'--r','MarkerIndices',1:50:Nx,'linewidth',3);
title('Lecture #25 Example','fontsize',16,'fontweight','bold');
xlabel('X [cm]','fontsize',14,'fontweight','bold');
ylabel('u(X,t)   [^{\circ}C]',...
    'fontsize',14,'fontweight','bold');
grid on
set(gca,'fontsize',12,'fontweight','bold');
legend('t = 0','t = 1','t = 10');

%% Simple time-dependent plot 
Tmax = 15; % s
NT = 150;
T = linspace(0,Tmax,NT);
figure(2)
for n = 1:NT
   plot(X,u(X,T(n)),'-b','linewidth',2)
   title_str = sprintf('Lecture #25 Example, t = %5.3g',T(n));
   title(title_str,'fontsize',16,'fontweight','bold');
   xlabel('X','fontsize',14,'fontweight','bold');
   ylabel('u(X,T)','fontsize',14,'fontweight','bold');
   grid on
   set(gca,'fontsize',12,'fontweight','bold');
   axis([0 L -0.2 4.5]);
   %pause(Tmax/(NT-1));
end

%% Save the time-dependent plot as a Movie (still more fun - optional)
FRAMES(NT) = struct('cdata',[],'colormap',[]);
figure(3)
for n = 1:NT
    plot(X,u(X,T(n)),'-b','linewidth',3);
    title_str = ...
        sprintf('Lecture #25 Example, t = %g ',T(n));
    title(title_str,'fontsize',16,'fontweight','bold');
    xlabel('X','fontsize',14,'fontweight','bold');
    ylabel('u(X,T)','fontsize',14,'fontweight','bold');
    grid on
    set(gca,'fontsize',12,'fontweight','bold');
    axis([0 L -0.2 4.5]);
    drawnow %<<< ensure graphics pipeline is complete/"flushed"
    FRAMES(n) = getframe(gcf);
end

%% play the movie
fig = figure(4);
movie(fig,FRAMES,10); %<< last argument is frames-per-second

%% Or Write to AVI (something to send to your Mom - optional)
v = VideoWriter('TransientHeat.avi');
open(v);
for n = 1:NT
   writeVideo(v,FRAMES(n)); 
end
close(v);
%% Plot the temperature vs time in a 2D plot using the surf function
[XX,TT] = meshgrid(X,T);
figure(5)
surf(XX,TT,u(XX,TT),'edgecolor','none'); 
title('Lecture 25 Suface Plot Example',...
    'fontsize',18,'fontweight','bold');
xlabel('X [cm]','fontsize',16,'fontweight','bold');
ylabel('T [s]','fontsize',16,'fontweight','bold');
zlabel('u(X,T)  [^{o}C]','fontsize',16,...
    'fontweight','bold');