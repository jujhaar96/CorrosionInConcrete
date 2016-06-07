function value = Approxf3(n, x )
%APPROXF Summary of this function goes here
% this function is the 2 order approximation of  2 variable function.
value = 0;
c = size(n);
s=c(1,1);
for i = 1:s
    v = Interp(n(:,1), n(i,1), x(1,1));
    w = Interp(n(:,2), n(i,2), x(1,2));
    value = value + v*w*fun(n(i,1),n(i,2));
end



end

