n=100000;
d=0.15;
gamma=1;
tau=0.4;
mu=0.2;
alpha=0.1;
kappa=3;
sigma=0.5;
k=kappa+sigma;
theta=0.8;
y0=0.1;
z0=0.1;
eps0=10;
T=100;
for i=1:100
    res(i,:)=CoEv_MC_nf(n,d,gamma,tau,mu,alpha,kappa,theta,T,y0,z0,eps0);
end
figure
plot(linspace(0,T,800),max(res))
hold on
plot(linspace(0,T,800),min(res))
plot(linspace(0,T,800),mean(res))


tspan = [0, T]; 

% Initial conditions
x0 = 0.1; 
y0 = 0.1;
z0 = 10; 

% Define the ODE system as a function
ode_system = @(t, y) [
 ((theta + d*(1-theta))*(d*min(y(1),1) + (1-d)*min(y(2),1))*min(y(1),1) + (1-d)*(1-theta)*(d*min(y(1),1) + (1-d)*min(y(2),1) + mu*y(3) + alpha)*min(y(2),1))*(1-min(y(1),1)) - (1 - d*min(y(1),1) - (1-d)*min(y(2),1) + k - sigma)*(1 - (theta + d*(1-theta))*min(y(1),1) - (1-d)*(1-theta)*min(y(2),1))*min(y(1),1);
    ((1 - d*(1-theta))*(d*min(y(1),1) + (1-d)*min(y(2),1) + mu*y(3) + alpha)*min(y(2),1) + d*(1-theta)*(d*min(y(1),1) + (1-d)*min(y(2),1))*min(y(1),1))*(1-min(y(2),1)) - (1 - d*min(y(1),1) - (1-d)*min(y(2),1) + k - sigma)*(1 - d*(1-theta)*min(y(1),1) - (1 - d*(1-theta))*min(y(2),1))*min(y(2),1);
    (gamma*(1 - d*min(y(1),1) - (1-d)*min(y(2),1)) - tau) * y(3)
];


% Solve the ODEs
[t, Y] = ode45(ode_system, tspan, [x0, y0, z0]);

figure
plot(d*Y(:, 1)+(1-d)*Y(:, 2),Y(:,3))

figure
plot(t,d*Y(:, 1)+(1-d)*Y(:, 2))

%figure
%plot(Y(:, 2),Y(:,3))