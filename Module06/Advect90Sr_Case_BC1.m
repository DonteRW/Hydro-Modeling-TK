% Purpose: This is a simple model to perform a simulation of the advection
% of 90^Sr in an aquifer with known conductivity and a constant regional
% hydraulic gradient.

clear all; close all;

dx = 1; % Spatial step in m
dt = 1; % Time step in yr

disp(['Courant number = ',num2str(dt/dx)]);

xmax = 4000; % Extent of domain in m
tmax = 20000; % Duration of the simulation in yr

Ks =   10; % Saturated hydraulic conductivity in cm/day
dhdx = 30;% Regional hydraulic gradient in ft/mi

% Convert to SI units of m/yr
Ks = Ks/100*365.25;
dhdx = dhdx*(1/5280); % A conversion from ft/mi to ft/ft, which is dimensionless (i.e., is also the slope in m/m);

% Compute the Darcy velocity
q = Ks*dhdx;

x = 0:dx:xmax; % Vector of our spatial "grid"
t = 0:dt:tmax; % Vector of simulation times

% Set initial and boundary conditions
% Initial condition: concentration is 0 Bq/g everywhere except at x = 0,
% where it is 1 Bq/g
c0 = zeros(length(x),1);
c0(1) = 1;

% Boundary condition: concentration is 1 Bq/g at x = 0 for all times
csource = zeros(length(t));
csource(1) = 1;

% Create a container to store simulated values of concentration
C = zeros(length(x),length(t));
C(:,1) = c0;

for i=2:length(t) % This is the time loop
    
    % Get the initial values
    Cprev = C(:,i-1);
    
    % Create a container for the values of Cnext (concentration at the next
    % time step)
    Cnext = zeros(length(x),1);
    
    for j=1:length(x) % This is the spatial loop
        if(j==1)
            Cnext(j) = csource(i); 
        else
            Cnext(j) = Cprev(j) + q*(dt/dx)*(Cprev(j-1) - Cprev(j));
        end
    end
    
    %
    
    % Store the value of the concentration profile at timestep i
    C(:,i) = Cnext;
    
    if(mod(t(i),500)==0)
        figure(1);
        plot(x,C(:,i),'b-'); hold on;
    end
    
end

% Plot some stuff...
figure(2);

subplot(411);
plot(t,C(x==500,:),'b-');
ylabel('()^{90}Sr conc. [Bq/g]');
ylim([0 0.05]);

subplot(412);
plot(t,C(x==1000,:),'b-');
ylabel('()^{90}Sr conc. [Bq/g]');
ylim([0 0.05]);

subplot(413);
plot(t,C(x==2000,:),'b-');
ylabel('()^{90}Sr conc. [Bq/g]');
ylim([0 0.05]);

subplot(414);
plot(t,C(x==4000,:),'b-');
ylabel('()^{90}Sr conc. [Bq/g]');
ylim([0 0.05]);
xlabel('Time [yr]');














