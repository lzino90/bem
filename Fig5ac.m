close all
clear all


% Define the parameters
theta=0.2;
dd=[0.01:0.01:0.26];
M=length(dd);
k = 2;         
mu = 1;      
sigma = 0;   
gamma = 1;   
tau = 0.8;   
alpha = 0;
L=100;

% Define the time span for the simulation
tspan = [0, 5000]; 

% Initial conditions
options = odeset('Mass',M,'RelTol',1e-4,'AbsTol',1e-6);

for i=1:L
    x0(i) = rand; 
    y0(i) = rand; 
    z0(i) = 5*rand; 
end


for i=1:M
    d=dd(i);
    ode_system = @(t, y) [
    ((theta + d*(1-theta))*(d*min(y(1),1) + (1-d)*min(y(2),1))*min(y(1),1) + (1-d)*(1-theta)*(d*min(y(1),1) + (1-d)*min(y(2),1) + mu*y(3) + alpha)*min(y(2),1))*(1-min(y(1),1)) - (1 - d*min(y(1),1) - (1-d)*min(y(2),1) + k - sigma)*(1 - (theta + d*(1-theta))*min(y(1),1) - (1-d)*(1-theta)*min(y(2),1))*min(y(1),1);
    ((1 - d*(1-theta))*(d*min(y(1),1) + (1-d)*min(y(2),1) + mu*y(3) + alpha)*min(y(2),1) + d*(1-theta)*(d*min(y(1),1) + (1-d)*min(y(2),1))*min(y(1),1))*(1-min(y(2),1)) - (1 - d*min(y(1),1) - (1-d)*min(y(2),1) + k - sigma)*(1 - d*(1-theta)*min(y(1),1) - (1 - d*(1-theta))*min(y(2),1))*min(y(2),1);
    (gamma*(1 - d*min(y(1),1) - (1-d)*min(y(2),1)) - tau) * y(3)
    ];
    count=0;
    for j=1:L
        [t, Y] = ode45(ode_system, tspan, [x0(j), y0(j), z0(j)]);
        z = Y(:, 3);
        if max(z)<=50
            mx_t(j)=max(z);
            av_t(j)=mean(z);
            count=count+1;
        end
    end
    display(strcat('Simulations:',num2str(round(i*100/M)),'%'))
    mx(i)=sum(mx_t)/count;
    av(i)=sum(av_t)/count;
    clear mx_t av_t
end
        
figure
plot(dd,mx)
figure
plot(dd,av)



% Define the parameters
theta=0.8;
dd=[0.01:0.01:0.26];
M=length(dd);
k = 2;         
mu = 1;      
sigma = 0;   
gamma = 1;   
tau = 0.8;   
alpha = 0;
L=3;

% Define the time span for the simulation
tspan = [0, 5000]; 

% Initial conditions
options = odeset('Mass',M,'RelTol',1e-4,'AbsTol',1e-6);

for i=1:L
    x0(i) = 0.1+0.8*rand; 
    y0(i) = 0.1+0.8*rand; 
    z0(i) = 5*rand; 
end


for i=1:M
    d=dd(i);
    ode_system = @(t, y) [
    ((theta + d*(1-theta))*(d*min(y(1),1) + (1-d)*min(y(2),1))*min(y(1),1) + (1-d)*(1-theta)*(d*min(y(1),1) + (1-d)*min(y(2),1) + mu*y(3) + alpha)*min(y(2),1))*(1-min(y(1),1)) - (1 - d*min(y(1),1) - (1-d)*min(y(2),1) + k - sigma)*(1 - (theta + d*(1-theta))*min(y(1),1) - (1-d)*(1-theta)*min(y(2),1))*min(y(1),1);
    ((1 - d*(1-theta))*(d*min(y(1),1) + (1-d)*min(y(2),1) + mu*y(3) + alpha)*min(y(2),1) + d*(1-theta)*(d*min(y(1),1) + (1-d)*min(y(2),1))*min(y(1),1))*(1-min(y(2),1)) - (1 - d*min(y(1),1) - (1-d)*min(y(2),1) + k - sigma)*(1 - d*(1-theta)*min(y(1),1) - (1 - d*(1-theta))*min(y(2),1))*min(y(2),1);
    (gamma*(1 - d*min(y(1),1) - (1-d)*min(y(2),1)) - tau) * y(3)
    ];
    count=0;
    for j=1:L
        [t, Y] = ode45(ode_system, tspan, [x0(j), y0(j), z0(j)]);
        z = Y(:, 3);
        if max(z)<=50
            mx_t(j)=max(z);
            av_t(j)=mean(z);
            count=count+1;
        end
    end
    display(strcat('Simulations:',num2str(round(i*100/M)),'%'))
    mx(i)=sum(mx_t)/count;
    av(i)=sum(av_t)/count;
    clear mx_t av_t
end
        
figure
plot(dd,mx)
figure
plot(dd,av)