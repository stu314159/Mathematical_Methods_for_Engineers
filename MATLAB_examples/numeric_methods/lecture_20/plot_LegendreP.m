% plot_LegendreP.m

% create an attractive plot of the Legendre Polynomials over the range
% [-1,1].

clear
clc
close all

P = 5; % < - P must be  0 > P <= 5 for this simple plot script

Pn = cell(P+1,1);
Pn{1} = @(x) 1;
Pn{2} = @(x) x;

for n = 2:P
    Pn{n+1} = @(x) (2*(n-1)+1)*x.*Pn{n}(x)./((n-1)+1) - (n-1)*Pn{n-1}(x)./((n-1)+1);
end

x = -1:.01:1;
colorspec = {'b','r','g','c','m','y'};

for n = 1:(P+1)
   
   plot(x,Pn{n}(x),colorspec{n},'linewidth',2);  
   axis([-1 1 -1.2 1.2])
   if n == 1
       hold on
   end
end

grid on
title(sprintf('First %i Legendre Polynomials.',P),...
    'fontsize',16,'fontweight','bold');
hold off