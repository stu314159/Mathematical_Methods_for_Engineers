% gauss_basis_test.m
%% Purpose

clear
clc
close all

f = @(x) exp(-(x.^2));
a = 0;
b = 3;
P = 4;

% get the gauss_quadrature basis
[xgl,wgl] = gauss_basis(P);

% variable change to scale interval to [-1,1] as required by Gaussian
% Integration
xT = @(t) ((b-a)*t + a + b)/2;
Jac = (b - a)/2;

% perform the numeric integration
intF = f(xT(xgl))*wgl*Jac;

% get result from MATLAB's built-in for comparison
% take this to be the "exact" value of the definite integral
intF_builtIn = integral(f,a,b);

fprintf('Result for my Gauss Quadrature for P = %i is: %18.16f \n',P,intF);
fprintf('Result for built in function integral = %18.16f \n',intF_builtIn);
fprintf('Relative difference = %g \n',abs(intF - intF_builtIn)/intF_builtIn);



