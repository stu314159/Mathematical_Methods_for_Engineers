%% Lecture 18 Demo
clear
clc
close 'all'

%%
% 
%  In this example, function and its derivative are represented with
%  Lagrange interpolants.  A choice is also provided between two sets of
%  sample points: a) evenly spaced points; and b) unevenly spaced points
%  placed at "Chebychev points".  The superiority of "Chebychev points" can
%  be seen by comparing the interpolation (and derivative) quality for,
%  especially, function 3 which is a scaled version of "the Witch of Agnesi" 
%  which is a function known to be interpolated poorly over evenly spaced
%  sample points due to the "Runge" phenomena
% 

%%

xMin = -2;
xMax = 2;
Nx = 51;

sample_set = 2;
% 1 = evenly spaced
% 2 = Chebychev points 

switch sample_set
    
    case 1
        X = linspace(xMin,xMax,Nx);
        
    case 2
        % map domain to [-1,1]
        xT = @(t) ((xMax - xMin)*t + xMin + xMax)/2;
        Tch = @(n) cos(((2*(1:n))-1)*pi./(2*n));
        
        xP = Tch(Nx);%<- get Nx Chebychev points
        X = xT(xP);%<-- map them to the domain of interest
end



fun_select = 3;
% 1 - x^3
% 2 - cos(kx)<-- k is variable
% 3 - "Witch of Agnesi"

switch fun_select
    
    case 1
        f = @(x) x.^3;
        df = @(x) 3.*(x.^2);
        
    case 2
        k = 5;
        f = @(x) cos(k*x);
        df = @(x) -k*sin(k*x);
        
    case 3
        f = @(x) 1./(1+25*x.^2);
        df = @(x) -50*x./((1+25*x.^2).^2);
        
end

Y = f(X);

dF_Lagrange = genLagrangeInterpDeriv(X,Y);

%% Plot the results
X_plt = linspace(xMin,xMax,1000);

figure(1)

subplot(2,1,1)
plot(X_plt,f(X_plt))
title('Lagrange Interpolation and Differentiation',...
    'FontSize',14,'FontWeight','bold');
%axis([xMin xMax -0.05 1.05]);

subplot(2,1,2)
plot(X_plt,dF_Lagrange(X_plt),'-b',...
    X,df(X),'sr');
axis([xMin xMax -2 2]);
xlabel('X','FontSize',12,'FontWeight','bold');


%% Local Functions
function dF = genLagrangeInterpDeriv(X,Y)
% function dF = genLagrangeInterpDeriv(X,Y) generates the 
% derivative of a function using Lagrange Polynomial 
% interpolation.
% Inputs
% X = x-values of the function
% Y = f(X) for some function
%
% Outputs
% dF - a function handle with the derivative of 
% the Lagrange interpolant

n = length(X);

dF = @(x) 0; %<-- start with zero, the additive identity

for i = 1:n
    dLi = @(x) 0;
    for k = 1:n
        if k ~= i
            dLp = @(x) 1;%<-- start with 1 (multiplicitive identity)
            for j = 1:n
                if ((j ~= i) && (j ~= k))
                    dLp = @(x) dLp(x).* ...
                        (x - X(j))./(X(i) - X(j));
                end
            end
            dLi = @(x) dLi(x) + ...
                (1./(X(i) - X(k))).*dLp(x);
        end
        
    end
    dF = @(x) dF(x) + dLi(x)*Y(i);
end
end