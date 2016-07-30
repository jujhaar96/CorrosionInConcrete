function X = Interp( n, j, x )
%INTERP Summary of this function goes here
%   Detailed explanation goes here
s= length(n);
X=1;
for i = 1:s
    if n(i) ~= j
    X= X*(x-n(i))/(j-n(i));
    end
end

