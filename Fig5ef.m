close all
clear all


% Define the parameters
th=[0.1:.02:.98];
dd=[0.02:0.02:0.26];
M=length(th);
R=length(dd);
d = 0.2;       
k = 2;         
mu = 1;      
sigma = 0;   
gamma = 1;   
tau = 0.8;   
alpha = 0;


% Define the time span for the simulation
tspan = [0, 5000]; 

% Initial conditions
options = odeset('Mass',M,'RelTol',1e-4,'AbsTol',1e-6);
L=100;
for i=1:L
    x0(i) = rand; 
    y0(i) = rand; 
    z0(i) = 5*rand; 
end


for r=1:R
    d=dd(r);
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
    display(strcat('Simulations:',num2str(round((i+(r-1)*M)*100/(M*R))),'%'))
    mx(i,r)=sum(mx_t)/count;
    av(i,r)=sum(av_t)/count;
    clear mx_t av_t
    end
end

%plot 4e
figure
surf(th,dd,mx')
view(0,90)

%plot 4f
figure
surf(th,dd,av')
view(0,90)