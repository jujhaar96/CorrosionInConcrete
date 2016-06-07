function Y = fun( x, z )
%FUN Summary of this function goes here
%   Detailed explanation goes here
si= max(size(x));
for i = 1:si

    Y(i,1) = x(i,1)^2 + z(i,1)^2;
end
end

