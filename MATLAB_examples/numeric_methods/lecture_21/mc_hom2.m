% mc_hom2.m
% The purpose of this script is to use the Monte Carlo method to integrate
% the function exp(-x^2) from [-3,3] and compare its performance with other
% numerical integration schemes.

% this uses the simple "hit-or-miss" algorithm discussed in section 2.1 of
% the "lecture 9" pdf in this directory.


clear
clc
close all

f = @(t) sin(1./(-t)).^2;

% this will be a convergence study, so I will need a range of sample sizes
N = 10:2:32;

xMin = 0;
xMax = 1;
yMin = 0;
yMax = 1;
trialArea = (xMax - xMin)*(yMax - yMin);

% estimate the exact value using built-in integrate function
%intExact = integral(f,xMin,xMax);
intExact = 0.6734568; % from Lecture 9

relErr = NaN(1,length(N));
trialInt = NaN(1,length(N));

ctr = 0;

% for *very* large simulations, break the problem into chuncks
chunkSize = 2^24;

tic;
for n = 2.^N
    ctr = ctr+1; % increment my counter
    
    % give an indication of progress
    fprintf('Estimating integral by Monte Carlo integration using %i samples.\n',n);
    
    numChunks = floor(n/chunkSize);
    numRemainder = n - (numChunks*chunkSize);
    
    intAccum = 0;
    for k = 1:numChunks
        % take the sample
        x = xMin + (xMax - xMin).*rand(chunkSize,1);
        y = yMin + (yMax - yMin).*rand(chunkSize,1);
        intAccum = intAccum + sum(y <=f(x));
    end
    
    if numRemainder > 0
        x = xMin + (xMax - xMin).*rand(numRemainder,1);
        y = yMin + (yMax - yMin).*rand(numRemainder,1);
        intAccum = intAccum + sum(y <= f(x));
    end
    
    trialInt(ctr) = (intAccum/n)*trialArea;
    relErr(ctr) = norm(trialInt(ctr) - intExact,2)/norm(intExact,2);
    
end
intTime = toc;
fprintf('Average of %g samples per second.\n',sum(2.^N)/intTime);
% plot the results

figure
loglog(2.^N,relErr,'-b','linewidth',2);
hold on
loglog(2.^N,sqrt(100)./(sqrt(2.^N)),'--r','linewidth',2);
title('Convergence of Monte Carlo Integration',...
    'fontsize',16,'fontweight','bold');
xlabel('Number of Sample Points','fontsize',14,'fontweight','bold');
ylabel('Relative Error','fontsize',14,'fontweight','bold');
legend('MC Integration','\theta(n^{-1/2}) Convergence',...
    'location','best');
grid on
set(gca,'fontsize',12,'fontweight','bold');