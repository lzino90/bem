function [y,z,eps,t]=CoEv_MC_ws(n,d,gamma,tau,mu,alpha,kappa,theta,T,y0,z0,eps0) %DA FARE
A=eye(n)+WS(n,20,.3);
A=ER(n,20/(n-1))+eye(n);
%A=ones(n);
for i=1:n
    nei{i}=find(A(i,:));
end
deg=sum(A,2);
nd=round(n*d);
na=n-nd;
for i=1:nd
    neis{i}=find(A(i,1:nd));
end
for i=nd+1:n
    neis{i}=find(A(i,nd+1:n));
end
for i=1:n
    degs(i)=length(neis{i});
end
y=round(y0*nd)/nd;
z=round(z0*na)/na;
eps=eps0;

x=zeros(n,1);
x(randperm(nd,y*nd))=1;
x(nd+randperm(na,y*na))=1;

t=1;
u1=zeros(n,1);
u0=zeros(n,1);
m=A*x./deg;

for i=1:nd
    u1(i)=m(i);
    u0(i)=1-m(i)+kappa;
end

for i=nd+1:n
    u1(i)=m(i)+mu*eps+alpha;
    u0(i)=1-m(i)+kappa;
end
time=0;

while time(t)<T 
    t=t+1;
    y(t)=y(t-1);
    z(t)=z(t-1);
    M=max(max(u0),max(u1));
    time(t)=time(t-1)+randexp(n*M);
    eps(t)=eps(t-1)*exp((gamma*(1-d*y(t-1)-(1-d)*z(t-1))-tau)*(time(t)-time(t-1)));
    i=randi(n);
    if rand<theta %same community
        j=neis{i}(randi(degs(i)));
    else %general
        j=nei{i}(randi(deg(i)));
    end
    if x(i)==0 && x(j)==1
        if rand<u1(j)/M
            x(i)=1;
            if i<=nd
                y(t)=y(t-1)+1/nd;
                z(t)=z(t-1);
            else
                y(t)=y(t-1);
                z(t)=z(t-1)+1/na;
            end
            m=A*x./deg;
            for i=1:nd
                u1(i)=m(i);
                u0(i)=1-m(i)+kappa;
            end
            for i=nd+1:n
                u1(i)=m(i)+mu*eps(t)+alpha;
                u0(i)=1-m(i)+kappa;
            end
        end
    elseif x(i)==1 && x(j)==0
        if rand<u0(j)/M
            x(i)=0;
            if i<=nd
                y(t)=y(t-1)-1/nd;
                z(t)=z(t-1);
            else
                y(t)=y(t-1);
                z(t)=z(t-1)-1/na;
            end
            m=A*x./deg;
            for i=1:nd
                u1(i)=m(i);
                u0(i)=1-m(i)+kappa;
            end
            for i=nd+1:n
                u1(i)=m(i)+mu*eps(t)+alpha;
                u0(i)=1-m(i)+kappa;
            end
        end
    end
end
            
 
[y,t]=reducev2(y,time,800);
[z,t]=reducev2(z,time,800);
[eps,t]=reducev2(eps,time,800);

  % figure
  %plot(y,eps,'r')
  %figure
   % plot(z,eps,'g')




    %hold on
  %  plot(t,z,'g')
  %  figure
 % plot(t,eps,'b')
end

