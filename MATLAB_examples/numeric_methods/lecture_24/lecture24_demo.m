%% Lecture 24 MATLAB Demo
clear
clc
close 'all'

%% Example #1 Xenon Transient

% Compute the time-history of Xenon concentration during a reactor power
% transient.  

% Given Data:

% core data
nominalFlux = 2e14; % n/cm^2-sec -- Average neutron flux
Sigma_F = .0452; %1/cm -- Macroscopic Fission Cross-Section in the core. 
%(based on typical U-235 atomic density atomic density and microscopic
%thermal neutron fission cross section for U-235)

% Iodine-135  and Tellurium-135 Nuclear data
gamma_Te = 0.061; % fission yield for Tellurium-135 (decays with 19-second half life to I-135)
lambda_I = 2.9173e-5; % 1/sec -- decay constant for I-135. Corresponds to a 6.6 hour half-life

% Xenon-135 Nuclear data
sigma_a_Xe = 2.6e6*10^(-24); % cm^2 -- Microscopic thermal neutron absorption cross section for Xe-135
gamma_Xe = 0.003; % fission yield for Xe-135 resulting from thermal neutron fission of U-235
lambda_Xe = 2.1185e-5; % 1/sec -- decay constant for Xe.  Corresponds to a 9.1 hour half-life


%% Time discretization
tStart = 0; % sec - time start
tEnd = 160*3600; %sec -  time end (160 hours)
numTs = 50000; % no fewer than 1 time step per hour.
tSpace = linspace(tStart,tEnd,numTs);
dT = tSpace(2)-tSpace(1);

%% Construct data arrays and provide initial conditions

P = nan(2,numTs+1); % an extra column for the initial data
P(1,1) = 0; % initial I-135 concentration;
P(2,1) = 0; % initial Xe-135 concentration;
pfrac = nan(1,numTs); % power fraction

%% Describe power history
% a start-up, run, then shutdown.
flux = @(t) nominalFlux*power_profile_Xe(t);

for ts = 1:numTs
    
    % say something comforting about program progress
    if mod(ts,10000)==0
        fprintf('Commencing time step %i.\n',ts);
    end
    
    % get the rate of change of poison concentrations
    % use a simple Forward-Euler time-stepping scheme
    dP = nan(2,1);
    T = tSpace(ts);
    pfrac(ts) = power_profile_Xe(T);
    % update Iodine concentration
    dP(1) = gamma_Te*Sigma_F*flux(T) - lambda_I*P(1,ts);  
    
    % update Xenon concentration
    dP(2) = gamma_Xe*Sigma_F*flux(T) + lambda_I*P(1,ts) ...
        - P(2,ts)*sigma_a_Xe*flux(T) - lambda_Xe*P(2,ts);
    
    % update the total poison concentration.
    P(:,ts+1) = P(:,ts) + dP*dT;
    
end

%% Plot your results
figure(1)
subplot(2,1,1)
plot(tSpace/3600,pfrac,'linewidth',2);
axis([0 160 -0.1 1.1]);
title('Power Profile')
set(gca,'fontsize',14,'fontweight','bold');
ylabel('Power Fraction',...
    'FontWeight','bold','FontSize',14);

subplot(2,1,2)
h = semilogy(tSpace/3600,P(2,1:(end-1)));
set(h,'linewidth',2);
set(gca,'fontsize',14,'fontweight','bold');
axis([0 160 2*10^14 2*10^16])
grid on
xlabel('Time (h)','FontWeight','bold','FontSize',14)
ylabel('Xenon-135 (at/cm^3)',...
    'FontWeight','bold','FontSize',14);
title('Xenon-135 Concentration','FontSize',14,'FontWeight','bold')

%% Example #2
x0 = 0; % initial displacement
v0 = 0; % initial velocity

tMax = 5; % seconds
nT = 100000;
T = linspace(0,tMax,nT);
w = nan(2,nT);
w(:,1) = [x0;v0]; % set initial conditions
dt = T(2)-T(1);
for t = 1:(nT-1)
    % say something comforting about program progress
    if mod(t,10000)==0
        fprintf('Commencing time step %i.\n',t);
    end
    dw = ex2(T(t),w(:,t));
    w(:,t+1) = w(:,t) + dw*dt;
end

x = w(1,:); %ft, displacement
v = w(2,:);

%% plot the results
figure(2)
plot(T,x,'-c','linewidth',2);
xlabel('Time [sec]','fontsize',14,'fontweight','bold');
ylabel('Displacement [ft]','fontsize',14,...
    'fontweight','bold');
title('Rocket Displacement vs. Time',...
    'fontsize',16,'fontweight','bold');
grid on
set(gca,'fontsize',10,'fontweight','bold');


%% Local Functions
function p = power_profile_Xe(t)
% function p = power_profile_Xe(t) - function, input a time t in seconds and
% returns a power level [0,1] giving the percent full power

tH = t/3600; % convert seconds to hours.
if tH < 60
    p = 1;
elseif tH>=60 && tH<80
    p = 0;
else
    p = 1;
end
end

function dw = ex2(t,w)
% define m(t)
mo = 100; % slugs, initial mass of fuel
m_dot = 1; % slug/s, fuel burn rate
m = @(t) mo - m_dot*t;% slugs, fuel mass

k = 2e5; % lb/ft, spring coefficent
c = 500; % slug/s, damping coefficient
T = 10000; % lb, thrust

% initialize and compute the output
dw = nan(2,1);
dw(1) = w(2);
dw(2) = T./m(t) - c*w(2)./m(t) - k*w(1)./m(t);
end
