%% Lecture 22 MATLAB Examples
clear
clc
close 'all'

%% trapz example #1

y = [0 0.3 0.5 1 1.5 2 2.5 3 4 5]; % m, depth
v = [0 0.4 0.5 0.56 0.6 0.63 0.66 0.68 0.71 0.74];% m/s, speed

h = 10; % m, width of channel.

Q = h*trapz(y,v); % m^3/s, volumetric flow rate
fprintf('Q = %g m^3/s \n',Q);

%% integral example #1
% definite single integral
f = @(x) x.*exp((x-2).^4);
a = 0; b = 4;
int_f = integral(f,a,b);
fprintf('I = %g \n',int_f);

%% integral example #2
% improper single integrals
f2 = @(x) 1./x.^2;

I2 = integral(f2,1,inf);
fprintf('Integral is: %g \n',I2);
I2_exact = 1;
true_error = abs(I2 - I2_exact);
true_rel_error = abs(I2 - I2_exact)/abs(I2_exact);
fprintf('Ex #2, True error: %g \n',true_error);
fprintf('Ex #2, True relative error: %g \n',true_rel_error);

%% Integral example #3
% improper integral with singularity at an end point

f3 = @(x) 1./sqrt(3-x);

I3 = integral(f3,0,3);
I3_exact = 2*sqrt(3);
fprintf('I3 = %g \n',I3);
true_rel_error = abs(I3 - I3_exact)/abs(I3_exact);
fprintf('True relative error: %g \n',true_rel_error);


%% Multiple Integral
fun = @(x,y,z) sqrt(x.^2+y.^2+z.^2).*exp(-(x.^2+y.^2+z.^2));
I4 = integral3(fun,-inf,inf,-inf,inf,-inf,inf);

I4_exact = 2*pi;

I4_rel_error = abs(I4 - I4_exact)/abs(I4_exact);
fprintf('I4 = %g \n',I4);
fprintf('I4 relative error: %g \n',I4_rel_error);

%% Multiple Integral example #2
R = 1; H = 10;
P = @(r,z) besselj(0,2.405*r/R).*cos(pi*z/H);

totP = 2*pi*integral2(@(r,z) P(r,z).*r,0,R,-H/2,H/2);
V = pi*(R.^2)*H;
AvgP = totP/V;
% peak power is just 1*C
fprintf('Peak/Avg power = %g \n',1./AvgP);



