function [ Y ] = Error( A, B )
%ERROR Summary of this function goes here
%   input approx and actual function
X=0;
Z=0;

si= max(size(A));
for i = 1:si
    X=X+(A(i,1)-B(i,1))^2;
end

for i = 1:si
    Z=Z+(B(i,1))^2;
end
Y=(X/Z)^0.5
end

