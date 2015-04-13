% BSBM_wrapper.m
clear all; close all;

OutputFile = 'FortRenoBSBM_BenchmarkRun.mat';

dt = 0.05;

% Load forcing data
load('SCAN2022_FortRenoOK.mat');

nSteps = length(P)/dt;
t = dt:dt:(length(P));

% Define parameters and store in a vector 'phi'
beta  = 2;
V1max = 8;
k12   = 0.045;
V2max = 20;
k23   = 0.035;
V3max = 150;
kb    = 0.001;

phi(1) = beta;
phi(2) = V1max;
phi(3) = k12;
phi(4) = V2max;
phi(5) = k23;
phi(6) = V3max;
phi(7) = kb;

% Create initial states
V10 = 3;
V20 = 10;
V30 = 75;

V0(1,1) = V10;
V0(2,1) = V20;
V0(3,1) = V30;

% Create container to store states and runoff
V     = nan*ones(3,nSteps);
V1    = nan*ones(1,nSteps);
V2    = nan*ones(1,nSteps);
V3    = nan*ones(1,nSteps);
qtro  = nan*ones(1,nSteps);
qb    = nan*ones(1,nSteps);
qdro  = nan*ones(1,nSteps);
mberr = nan*ones(1,nSteps);

for i=1:nSteps
    
    % On the first iteration set the current condition to the initial
    % condition
    if(i==1)
        Vn = V0;
    else
        Vn = V(:,i-1);
    end
    
    % Set forcings
    day = ceil(i/(1/dt));
    
    U(1) = P(day);
    U(2) = 0.0;
    U(3) = 0.0;
    
    % Call the RK4 scheme
    [V(:,i),qdro(i),qb(i),qtro(i),mberr(i)] = rk4_bsbm(Vn,phi,U,dt);
    V1(i) = V(1,i);
    V2(i) = V(2,i);
    V3(i) = V(3,i);
    
end

% Make plots
figure(1);
subplot(311);
plot(t,V1);
ylabel('V_1 [mm]');
subplot(312);
plot(t,V2);
ylabel('V_2 [mm]');
subplot(313);
plot(t,V3);
ylabel('V_3 [mm]');
xlabel('Time [days]');

figure(2);
subplot(211);
plot(t,qdro,'r'); hold on;
plot(t,qb,'b');
plot(t,qtro,'k');
legend('Direct runoff','Base flow','Total runoff');
ylabel('Runoff rate [mm/d]');
subplot(212);
plot(t,mberr);
ylabel('Mass balance errors [mm]');

save(OutputFile,'V','V1','V2','V3','qdro','qb','qtro','mberr');
