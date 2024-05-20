close all
clear all


% Define the parameters
aa=[0:0.1:2];
mumu=[1:0.2:4];
M=length(aa);
R=length(mumu);
d = 0.2;
%k = 2;         
%mu = 1;      
sigma = 0;   
gamma = 1;   
tau = 0.8;   
alpha = 0;


% Define the time span for the simulation
tspan = [0, 200]; 

% Initial conditions
options = odeset('RelTol',1e-6,'AbsTol',1e-8,'NonNegative',1);
L=100;
for i=1:L
    y0(i) = rand; 
    z0(i) = rand*0.2; 
end


for r=1:R
    mu=mumu(r);
    for i=1:M
        k=4-aa(i);
        ode_system = @(t, y) [y(1)*(1-y(1))*(2*(1-d)*y(1)+mu*y(2)+alpha-k-1);
         (gamma*(1 - (1-d)*max(y(1),0)) - tau) * max(y(2),0)];
        count=0;
        for j=1:L
            [t, Y] = ode45(ode_system, tspan, [y0(j), z0(j)],options);
            z = Y(:, 2);
            if max(z)<=200
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
surf(aa,mumu,mx')
view(0,90)

%plot 4f
figure
surf(aa,mumu,av')
view(0,90)