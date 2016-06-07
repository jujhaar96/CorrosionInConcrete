function X = Interp3( n, j, x )
%INTERP Summary of this function goes here
%   Detailed explanation goes here
s= max(size(n));
X=1;
for i = 1:s
    if n(i,2) ~= j
    X= X*(x-n(i,2))/(j-n(i,2));
    end
end
