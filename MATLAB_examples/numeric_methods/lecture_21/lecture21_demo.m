%% EM486A Lecture 21 Monte Carlo Integration Demo

clear
clc
close 'all'

%% Choose function
fun_opt = 2;

switch fun_opt
    case 1
        f = @(t) exp(-t.^2);
    case 2
        f = @(t) sin(1./(-t)).^2;
        %intExact ~ 2.48308
    otherwise
        error('Invalid function choice!\n');
end

%% Establish the Bounding Box
xMin = -3;
xMax = 3;
yMin = 0;
yMax = 1;
trialArea = (xMax - xMin)*(yMax - yMin);

%% estimate the exact value using built-in integrate function
intExact = integral(f,xMin,xMax);

figure(1)
fplot(f,[xMin,xMax],'linewidth',3)
grid on
title('Function to be Integrated','fontsize',16,'fontweight','bold');
xlabel('X','fontsize',14,'fontweight','bold');
ylabel('F(X)','fontsize',14,'fontweight','bold');

%% Determine number of samples to take

n = 2^24;

% for *very* large simulations, break the problem into chunks
chunkSize = 2^24;

numChunks = floor(n/chunkSize);
numRemainder = n - (numChunks*chunkSize);


%% Hit-or-Miss Algorithm
fprintf('Hit or Miss Algorithm \n\n');
intAccum = 0;
tic;
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
int_time = toc;
MC_Int = (intAccum/n)*trialArea;
relErr = norm(MC_Int - intExact,2)/norm(intExact,2);

fprintf('(Hit or Miss) Estimated Integral: %15.14f \n',MC_Int);
fprintf('Relative Error: %g \n',relErr);
fprintf('(Hit or Miss) Avg number of samples per second: %g \n',n/int_time);


%% Mean Value Algorithm
fprintf('\nMean Value Algorithm \n');
intAccum = 0;
tic;
for k = 1:numChunks
    % take the sample
    x = xMin + (xMax - xMin).*rand(chunkSize,1);
    intAccum = intAccum + sum(f(x));
end

if numRemainder > 0
    x = xMin + (xMax - xMin).*rand(numRemainder,1);
    intAccum = intAccum + sum(f(x));
end
int_time = toc;
MC_Int = (intAccum/n)*trialArea;
relErr = norm(MC_Int - intExact,2)/norm(intExact,2);

fprintf('(Mean Value) Estimated Integral: %15.14f \n',MC_Int);
fprintf('Relative Error: %g \n',relErr);
fprintf('(Mean Value) Avg number of samples per second: %g \n',n/int_time);