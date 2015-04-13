clear all; close all;

x = 0:1:100; % Define the spatial node locations
t = 0:1:600; % Define the temporal vector

K = 3; % Saturated hydraulic conductivity

deltax = x(2) - x(1); % Back calculate delta x
deltat = t(2) - t(1); % Back calculate delta t
alpha = K*deltat/deltax^2; % Compute the coefficient alpha

Nx = length(x); % Get the length of the x-vector
Nt = length(t); % Get the length of the t-vector

Delta2 = diag(-2*ones(Nx,1)) + diag(ones(Nx-1,1),1) + diag(ones(Nx-1,1),-1); % Define the second differencing matrix
I = eye(Nx); % Get the Nx-dimensional identity matrix

hi = 8; % The initial pressure head (everywhere)
hf = 4; % The final pressure head (on the left boundary)

h = zeros(Nx,Nt); % Preallocate space for the pressure heads in the simulation

% Specify initial conditions
h0 = hi*ones(Nx,1);
h(:,1) = h0;

% Specify boundary conditions on the left
h_left = hi - (hi - hf)/(t(Nt)-t(1))*t;

% Plot the initial conditions
figure(1);
plot(x,h(:,1),'b-'); hold on;

% Construct the A matrix
A = I - alpha*Delta2;
A(1,:) = [1,zeros(1,(Nx-1))];    
A(Nx,:) = [zeros(1,(Nx-1)),1];

for i=2:Nt
    
    % Get the initial head at this time step
    hinit     = h(:,i-1);
    
    % Solve for hnext
    hnext = pinv(A)*hinit;
    
    hnext(1) = h_left(i);
    hnext(Nx) = hi;
    
    % Store hnext in the matrix h
    h(:,i) = hnext;
    
    if(mod(i,20)==0) % plot every 20th head distribution
        plot(x,h(:,i),'b-'); hold on;
        xlabel('Distance');
        ylabel('Head');
    end
    
end

figure(2);
surf(x,t,h'); shading interp;
xlabel('Distance');
ylabel('Time');
zlabel('Pressure head');
