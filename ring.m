function A = ring(n,m)
%A = ring(n,m) generates a ring with network size n and m number of
%connections (per side)

A=zeros(n); 

for i=1:n-m
    for j=i+1:i+m
        A(i,j)=1;
        A(j,i)=1;
    end
end
for i=n-m+1:n
     for j=i+1:n
         A(i,j)=1;
         A(j,i)=1;
     end
     for j=1:m-n+i
         A(i,j)=1;
         A(j,i)=1;
     end
end
%More complactly
% for i=1:m
%     A=A+diag(ones(n-i,1),i)+diag(ones(n-i,1),-i)+diag(ones(i,1),n-i)+diag(ones(i,1),-n+i); %we create the diagonal bands
% end



end