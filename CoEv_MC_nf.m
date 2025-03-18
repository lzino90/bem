function res=CoEv_MC_nf(n,d,gamma,tau,mu,alpha,kappa,theta,T,y0,z0,eps0)
nd=round(n*d);
na=n-nd;
y=round(y0*nd)/nd;
z=round(z0*na)/na;
eps=eps0;

t=1;

r01d=(theta+d*(1-theta))*(d*y(t)+(1-d)*z(t))*y(t) + (1-d)*(1-theta)*(d*y(t)+(1-d)*z(t)+mu*eps(t)+alpha)*z(t);
r01a=(1-d*(1-theta))*(d*y(t)+(1-d)*z(t)+mu*eps(t)+alpha)*z(t)+d*(1-theta)*(d*y(t)+(1-d)*z(t))*y(t);
r10d=(1-d*y(t)-(1-d)*z(t)+kappa)*(1-(theta+d*(1-theta))*y(t)-(1-d)*(1-theta)*z(t));
r10a=(1-d*y(t)-(1-d)*z(t)+kappa)*(1-d*(1-theta)*y(t)-(1-d*(1-theta))*z(t));



time=0;
energy=n*(d*(1-y(t))*r01d+d*y(t)*r10d+(1-d)*(1-z(t))*r01a+(1-d)*z(t)*r10a);

while time(t)<T && energy>0
    t=t+1;
    y(t)=y(t-1);
    z(t)=z(t-1);
    time(t)=time(t-1)+randexp(energy);
    eps(t)=eps(t-1)*exp((gamma*(1-d*y(t-1)-(1-d)*z(t-1))-tau)*(time(t)-time(t-1)));
    prob=n*[d*(1-y(t))*r01d,d*y(t)*r10d,(1-d)*(1-z(t))*r01a,(1-d)*z(t)*r10a]/energy;
    temp=randvett(sommacumulativa(prob));
    switch temp
        case 1
            y(t)=y(t-1)+1/nd;
        case 2
            y(t)=y(t-1)-1/nd;
        case 3
            z(t)=z(t-1)+1/na;
        otherwise
            z(t)=z(t-1)-1/na;
    end
r01d=max(0,(theta+d*(1-theta))*(d*y(t)+(1-d)*z(t))*y(t) + (1-d)*(1-theta)*(d*y(t)+(1-d)*z(t)+mu*eps(t)+alpha)*z(t));
r01a=max(0,(1-d*(1-theta))*(d*y(t)+(1-d)*z(t)+mu*eps(t)+alpha)*z(t)+d*(1-theta)*(d*y(t)+(1-d)*z(t))*y(t));
r10d=max(0,(1-d*y(t)-(1-d)*z(t)+kappa)*(1-(theta+d*(1-theta))*y(t)-(1-d)*(1-theta)*z(t)));
r10a=max(0,(1-d*y(t)-(1-d)*z(t)+kappa)*(1-d*(1-theta)*y(t)-(1-d*(1-theta))*z(t)));
energy=max(0,n*(d*(1-y(t))*r01d+d*y(t)*r10d+(1-d)*(1-z(t))*r01a+(1-d)*z(t)*r10a));
end
            
 
 [y,t]=reducev2(y,time,800);
  [z,t]=reducev2(z,time,800);
   [eps,t]=reducev2(eps,time,800);

  % figure
  %plot(y,eps,'r')
  %figure
   % plot(z,eps,'g')

%  figure
 %   plot(d*y+(1-d)*z,eps,'g')  
    
  %  figure
   %     plot(t,d*y+(1-d)*z)  
    res=d*y+(1-d)*z;

    %hold on
  %  plot(t,z,'g')
  %  figure
 % plot(t,eps,'b')
end

