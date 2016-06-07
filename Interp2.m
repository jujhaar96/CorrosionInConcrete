function Y = Interp2( n, j, x)
%UNTITLED Summary of this function goes here
%   j and x should be of the same dimension. we haven't provided the check
%   here. LATER.
Y=1;
s= max(size(j));

for i = 1:s
Y = Y * Interp(n, j(i,1), x(i,1));
end
end

