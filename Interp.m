function X = Interp( n, j, x )
%INTERP Summary of this function goes here
%   Detailed explanation goes here
s= max(size(n));
X=1;
for i = 1:s
    if n(i,1) ~= j
    X= X*(x-n(i,1))/(j-n(i,1));
    end
end

