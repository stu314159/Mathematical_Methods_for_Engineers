% mc_hit_or_miss_schematic.m

clear
clc
close all

f = @(t) exp(-t.^2);

N = 1000;

xMin = -3;
xMax = 3;
yMin = 0;
yMax = 1;

trialArea = (xMax - xMin)*(yMax - yMin);

xt = xMin + (xMax - xMin).*rand(N,1);
yt = yMin + (yMax - yMin).*rand(N,1);

x = xMin:.01:xMax;

figure
plot(x,f(x),'linewidth',2);
hold on
plot(xt,yt,'.r','markersize',10);
set(gca,'fontsize',12,'fontweight','bold');
title('Hit or Miss Approach','fontsize',16,'fontweight','bold');
xlabel('x','fontsize',14,'fontweight','bold');
ylabel('y','fontsize',14,'fontweight','bold');
grid on
