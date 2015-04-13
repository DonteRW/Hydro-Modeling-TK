clear all; close all;

ti = 0;   % Initial time
tf = 600; % Final time

K = 3; % Saturated hydraulic conductivity

x = 0:1:100; % Define the spatial node locations
deltax = x(2) - x(1); % Back calculate delta x

alpha = 0.4; % DEFINE the coefficient alpha for stability
deltat = alpha*deltax^2/K; % Back calculate delta t based on alpha, deltax, and K 
t = ti:deltat:tf; % Define the t-vector

Nx = length(x); % Get the length of the x-vector
Nt = length(t); % Get the length of the t-vector

disp(['Number of time steps = ',int2str(Nt)]); % Display the number of time steps to contrast with implicit methods

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

for i=2:Nt
    
    % Get the initial head at this time step
    hinit     = h(:,i-1);
    
    % Solve for hnext
    hnext = hinit + alpha*Delta2*hinit;
    
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

% Make a surface plot of the simulated pressure heads
figure(2);
surf(x,t,h'); shading interp;
xlabel('Distance');
ylabel('Time');
zlabel('Pressure head');
