function y = modn(x,n)
%modn(x,n) computes x modulus n (the reminder of the division x:n), but it returns n instead of z is the
%reminder is 0
y=mod(x,n);
if y==0
    y=n;
end
end