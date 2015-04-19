% Purpose: This is a simple model to perform a simulation of the advection
% of 90^Sr in an aquifer with known conductivity and a constant regional
% hydraulic gradient.

clear all; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dx = 5; % Spatial step in m
D = 1; % Dispersivity of Strontium-90 in aquifer

Ks =   5; % Saturated hydraulic conductivity in cm/day
dhdx = 30; % Regional hydraulic gradient in ft/mi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xmax = 4000; % Extent of domain in m
tmax = 20000; % Duration of the simulation in yr

% Convert to SI units of m/yr
Ks = Ks/100*365.25;
dhdx = dhdx*(1/5280); % A conversion from ft/mi to ft/ft, which is dimensionless (i.e., is also the slope in m/m);

% Compute the Darcy velocity
q = Ks*dhdx;

% Back-calculate dt to maintain stability of the numerical solution
dt = dx;

flag = 0;
count = 1;
while (flag==0)

    % Check for stability and adjust if need be
    c = q*dt/dx; % Courant condition
    alpha = D*dt/dx^2; % Diffusion number

    if((c^2 <= 2*alpha)&&((alpha + c/4) <= (1/2)))
        flag = 1;
    else
        dt = dt/2; % Halve the calculated time step
    end
    
    count = count + 1;
end

% Output Courant condition
disp(['Courant condition = ',num2str(c)]);

% Output diffusion number
disp(['Diffusion number = ',num2str(D*dt/dx^2)]);

% Output time step
disp(['The time-step calculated for stability is = ',num2str(dt)]);


x = 0:dx:xmax; % Vector of our spatial "grid"
t = 0:dt:tmax; % Vector of simulation times

% Set initial and boundary conditions
% Initial condition: concentration is 0 Bq/g everywhere except at x = 0,
% where it is 1 Bq/g
c0 = zeros(length(x),1);
c0(1) = 1;

% Boundary condition: concentration is 1 Bq/g at x = 0 for all times
csource = zeros(length(t),1);
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
        % Finite difference solution to the 1-D ADE using the QUICK scheme
        if(j==1)
            Cnext(j) = csource(i);
        elseif(j==2) % At node 1, assume that the behavior at node -1 (i.e., j-2) is the same as node 1 at the previous time step
            Cnext(j) = Cprev(j) - c*((1/8)*Cprev(j) - (7/8)*Cprev(j-1) + (3/8)*Cprev(j) + (3/8)*Cprev(j+1)) ...
                + alpha*(Cprev(j-1) - 2*Cprev(j) + Cprev(j+1));            
        elseif(j==length(x))
            Cnext(j) = Cprev(j) - c*((1/8)*Cprev(j-2) - (7/8)*Cprev(j-1) + (3/8)*Cprev(j) + (3/8)*Cprev(j)) ...
                + alpha*(Cprev(j-1) - 2*Cprev(j) + Cprev(j));
        else 
            Cnext(j) = Cprev(j) - c*((1/8)*Cprev(j-2) - (7/8)*Cprev(j-1) + (3/8)*Cprev(j) + (3/8)*Cprev(j+1)) ...
                + alpha*(Cprev(j-1) - 2*Cprev(j) + Cprev(j+1));
        end
    end
    
    % Store the value of the concentration profile at timestep i
    C(:,i) = Cnext;
    
    if(mod(i,floor(500/dt))==0)
        figure(1);
        plot(x,C(:,i),'b-'); hold on;
    end
    
end

% Plot some stuff...
figure(2);

subplot(411);
plot(t,C(x==500,:),'b-');
ylabel('()^{90}Sr conc. [Bq/g]');
ylim([0 0.01]);

subplot(412);
plot(t,C(x==1000,:),'b-');
ylabel('()^{90}Sr conc. [Bq/g]');
ylim([0 0.01]);

subplot(413);
plot(t,C(x==2000,:),'b-');
ylabel('()^{90}Sr conc. [Bq/g]');
ylim([0 0.01]);

subplot(414);
plot(t,C(x==4000,:),'b-');
ylabel('()^{90}Sr conc. [Bq/g]');
ylim([0 0.01]);
xlabel('Time [yr]');














