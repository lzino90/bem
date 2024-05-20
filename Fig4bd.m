close all
clear all


% Define the parameters
th=[0.1:.01:.98];
d=0.1;
M=length(th);
k = 2;         
mu = 1;      
sigma = 0;   
gamma = 1;   
tau = 0.8;   
alpha = 0;
L=100;

% Define the time span for the simulation
tspan = [0, 10000]; 

% Initial conditions
options = odeset('RelTol',1e-4,'AbsTol',1e-6);

for i=1:L
    x0(i) = rand; 
    y0(i) = rand; 
    z0(i) = 5*rand; 
end


for i=1:M
    theta=th(i);
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
    display(strcat('Simulations:',num2str(round(i*50/M)),'%'))
    mx(i)=sum(mx_t)/count;
    mx_sd(i)=sqrt(sum((mx_t-mx(i)).^2)/count);
    av(i)=sum(av_t)/count;
    av_sd(i)=sqrt(sum((av_t-av(i)).^2)/count);
    clear mx_t av_t
end
        
figure
%plot(th,mx)
errorbar(th,mx,1.96*mx_sd)
figure
%plot(th,av)
errorbar(th,av,1.96*av_sd)



% Define the parameters
th=[0.1:.01:.98];
d=0.2;
M=length(th);
k = 2;         
mu = 1;      
sigma = 0;   
gamma = 1;   
tau = 0.8;   
alpha = 0;
L=100;

% Define the time span for the simulation
tspan = [0, 10000]; 

% Initial conditions
options = odeset('RelTol',1e-4,'AbsTol',1e-6);

for i=1:L
    x0(i) = rand; 
    y0(i) = rand; 
    z0(i) = 5*rand; 
end


for i=1:M
    theta=th(i);
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
    display(strcat('Simulations:',num2str(round(50+i*50/M)),'%'))
  mx(i)=sum(mx_t)/count;
    mx_sd(i)=sqrt(sum((mx_t-mx(i)).^2)/count);
    av(i)=sum(av_t)/count;
    av_sd(i)=sqrt(sum((av_t-av(i)).^2)/count);
    clear mx_t av_t
end
        
figure
%plot(th,mx)
errorbar(th,mx,1.96*mx_sd)
figure
%plot(th,av)
errorbar(th,av,1.96*av_sd)

