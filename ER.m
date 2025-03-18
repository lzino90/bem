function A = ER(n,p)
%A = ER(n,p) generates an Erdos–Renyi random graph with network size n and 
% link probability p

% When you write functions, it is often a good idea to put some checks,
% just to verify that inputs are correctly defined (e.g., that p is between
% 0 and 1, n positive integer,...)

A=zeros(n); 

%% Remark: How to simulate an event (in our case, generate a link) that occurs with a certain probability p? Generate a random number and if the random number is less than or equal to p, than the event occurs

for i=1:n
    for j=i+1:n
        if rand<=p
            A(i,j)=1;
            A(j,i)=1;
        end
    end
end

%% We can do it more efficiently as: A=(rand(n)<=p); A=triu(A,1)+triu(A,1)'

end