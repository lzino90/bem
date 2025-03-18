function A = WS(n,d,p)
m=d/2;
A=ring(n,m); 
for i=1:n
    for j=i+1:i+m 
        %I create an additional function to take care for links of the form (n,1)
        if A(i,modn(j,n))==1
            %check if we rewire it
            if rand<=p
                % rewire
                flag=0;
                while flag==0
                    S=randi(n);
                    if A(i,S)==0 && S~=i  %arc not present and not self loop
                        A(i,modn(j,n))=0; %cancel existing arc
                        A(modn(j,n),i)=0; %in both directions
                        A(i,S)=1; %create new arc
                        A(S,i)=1;
                        flag=1; %exit the loop
                    end  
                end
                %no instructions in the "else" are needed, since we keep the original arc in A
            end
        end
    end
end



end