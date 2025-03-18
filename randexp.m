function [x] = randexp(lambda,n)
if nargin==2
    lambda=lambda*ones(1,n);
end

if lambda<=0
    error('Non ha senso un parametro negativo');
end
h=length(lambda);
u=rand(1,h);
x=-log(u)./lambda;
end

