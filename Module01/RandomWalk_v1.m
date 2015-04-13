% Description: This is a one-dimensional random walk code
% Coder: Lejo Flores
% Date: 2015-01-22
% Assumptions:
%   - Starting position is 0
%   - Step lengths are normally distributed with a mean xbar and a standard
%     deviation sx

% Control variables:
Nwalkers = 100;
Nsteps = 1000;
xbar = 2;
sx = 5;
x0 = 0;

% Create a container to store walker location
x = zeros(Nsteps,Nwalkers);

for i=1:Nwalkers
    
    for j=1:Nsteps
        if(j==1)
            x(j,i) = x0;
        else
            x_step = xbar + sx*randn(1);
            x(j,i) = x(j-1,i) + x_step;
        end
    
    end
    plot(x(:,i)); hold on;
    
end

xlabel('Step');
ylabel('Position');
