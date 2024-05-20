close all
clear all

purple = [117 26 255]./255;


% Define the parameters
theta = 0.8;   
d = 0.1;       
k = 4;         
mu = 0.2;      
sigma = 0.5;   
gamma = 1;   
tau = 0.4;   
alpha = 0.1;


% Define the time span for the simulation
tspan = [0, 2000]; 

% Initial conditions
x0 = 0.001; 
y0 = 0.05;
z0 = 10; 

% Define the ODE system as a function
ode_system = @(t, y) [
 ((theta + d*(1-theta))*(d*min(y(1),1) + (1-d)*min(y(2),1))*min(y(1),1) + (1-d)*(1-theta)*(d*min(y(1),1) + (1-d)*min(y(2),1) + mu*y(3) + alpha)*min(y(2),1))*(1-min(y(1),1)) - (1 - d*min(y(1),1) - (1-d)*min(y(2),1) + k - sigma)*(1 - (theta + d*(1-theta))*min(y(1),1) - (1-d)*(1-theta)*min(y(2),1))*min(y(1),1);
    ((1 - d*(1-theta))*(d*min(y(1),1) + (1-d)*min(y(2),1) + mu*y(3) + alpha)*min(y(2),1) + d*(1-theta)*(d*min(y(1),1) + (1-d)*min(y(2),1))*min(y(1),1))*(1-min(y(2),1)) - (1 - d*min(y(1),1) - (1-d)*min(y(2),1) + k - sigma)*(1 - d*(1-theta)*min(y(1),1) - (1 - d*(1-theta))*min(y(2),1))*min(y(2),1);
    (gamma*(1 - d*min(y(1),1) - (1-d)*min(y(2),1)) - tau) * y(3)
];

% ode_system = @(t, y) [
%  ((theta + d*(1-theta))*(d*y(1) + (1-d)*y(2))*y(1) + (1-d)*(1-theta)*(d*y(1) + (1-d)*y(2) + mu*y(3) + alpha)*y(2))*(1-y(1)) - (1 - d*y(1) - (1-d)*y(2) + k - sigma)*(1 - (theta + d*(1-theta))*y(1) - (1-d)*(1-theta)*y(2))*y(1);
%     ((1 - d*(1-theta))*(d*y(1) + (1-d)*y(2) + mu*y(3) + alpha)*y(2) + d*(1-theta)*(d*y(1) + (1-d)*y(2))*y(1))*(1-y(2)) - (1 - d*y(1) - (1-d)*y(2) + k - sigma)*(1 - d*(1-theta)*y(1) - (1 - d*(1-theta))*y(2))*y(2);
%     (gamma*(1 - d*y(1) - (1-d)*y(2)) - tau) * y(3)
% ];

% Solve the ODEs
[t, Y] = ode45(ode_system, tspan, [x0, y0, z0]);

% Extract the solution components
 p = Y(:, 1);
 q = Y(:, 2);
 z = Y(:, 3);
 x = d*p +(1-d)*q;

% Plot the expression d*x_D + (1-d)*x_C on the x-axis and epsilon on the y-axis

 figure;
 plot(p, z, 'LineWidth', 1.5,'Color','r');
   xlim([0 1])
   ylim([0 60])
 xlabel('Responsible behaviour of deniers $x_{\mathcal{D}}$','interpreter', 'latex','FontSize', 23);
 ylabel('Environmental impact $\varepsilon$', 'interpreter', 'latex','FontSize', 23);
%title('Trajectory');
 grid on;

 figure;
 plot(q, z, 'LineWidth', 1.5,'Color','b');
   xlim([0 1])
   ylim([0 60])
 xlabel('Responsible behaviour of acknowledgers $x_{\mathcal{C}}$','interpreter', 'latex','FontSize', 23);
 ylabel('Environmental impact $\varepsilon$', 'interpreter', 'latex','FontSize', 23);
%title('Trajectory');
 grid on;

figure;
 plot(x, z, 'LineWidth', 1.5,'Color',purple);
   xlim([0 1])
   ylim([0 60])
 xlabel('Responsible behaviour of total population $x$','interpreter', 'latex','FontSize', 23);
 ylabel('Environmental impact $\varepsilon$', 'interpreter', 'latex','FontSize', 23);
%title('Trajectory');
 grid on;

% Create a 3D trajectory plot
 figure;
 plot3(q, p, z, 'LineWidth', 1.5,'Color','b');
  xlim([0 1])
  ylim([0 1])
 zlim([0 60])
xlabel('$$x_{\mathcal{C}}$$', 'interpreter', 'latex','FontSize', 20);
 ylabel('$x_{\mathcal{D}}$', 'interpreter', 'latex','FontSize', 20);
 zlabel('Environmental impact $$\varepsilon$$', 'interpreter', 'latex','FontSize', 20);
% title('3D Trajectory Plot');
 grid on;