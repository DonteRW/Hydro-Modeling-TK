clear all; close all;

x = 0:1:100; % Define spatial node locations
t = 0:1:600; % Define time vector

K = 3; % Saturated hydraulic conductivity

deltax = x(2) - x(1); % Back calculate deltax
deltat = t(2) - t(1); % Back calculate deltat
alpha = K*deltat/(2*deltax^2); % Calculate the coefficient alpha

Nx = length(x); % Get the number of spatial nodes
Nt = length(t); % Get the number of time steps

hi = 8; % Initial pressure head (everywhere)
hf = 4; % Final pressure head (on the left)

h = zeros(Nx,Nt); % Preallocate space to store the simulated pressure heads

% Specify initial conditions
h0 = 8*ones(Nx,1);
h(:,1) = h0;

% Specify boundary conditions
h_left = hi - (hi - hf)/(t(Nt)-t(1))*t;

% Plot the initial conditions
figure(1);
plot(x,h(:,1),'b-'); hold on;

% Construct the A matrix
A = diag((1+2*alpha)*ones(Nx,1)) + diag(-alpha*ones((Nx-1),1),-1) + diag(-alpha*ones((Nx-1),1),1);
A(1,:) = [1,zeros(1,(Nx-1))];
A(Nx,:) = [zeros(1,(Nx-1)),1];

% Construct the B matrix
B = diag((1-2*alpha)*ones(Nx,1)) + diag(alpha*ones((Nx-1),1),-1) + diag(alpha*ones((Nx-1),1),1);
B(1,:) = [1,zeros(1,(Nx-1))];
B(Nx,:) = [zeros(1,(Nx-1)),1];

for i=2:Nt
    
    % Get the initial head at this time step
    hinit = h(:,i-1);
    
    % Solve for hnext
    hnext = pinv(A)*B*hinit;
    hnext(1)  = h_left(i); 
    hnext(Nx) = hi;
    
    % Store hnext in the matrix h
    h(:,i) = hnext;
    
    if(mod(i,20)==0) % plot every 10th head distribution
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
