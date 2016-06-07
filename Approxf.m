function value = Approxf(n, x )
%APPROXF Summary of this function goes here
% this function is the 1 order approximation of  2 variable function.
value = 0;

a = unique(n(:,1));
b = unique(n(:,2));
c = length(a);
for i = 1:c
    v = Interp(a(:,1), a(i,1), x(1,1))*fun(a(i,1),0.5);
    w = Interp(b(:,1), b(i,1), x(1,2))*fun(0.5,b(i,1));
    value = value + v + w;
end

value = value - fun(0.5,0.5);


end

