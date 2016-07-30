function Y = fun( x )
%FUN Summary of this function goes here
%   Detailed explanation goes here
s= length(x);
Y=0;
for i = 1:(s-1)

    Y = Y+(x(i)*x(i+1));
end