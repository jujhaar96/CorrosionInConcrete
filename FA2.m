function Y = FA2( n, x)
%UNTITLED Summary of this function goes here
%   j and x should be of the same dimension. we haven't provided the check
%   here. LATER.
Y=zeros(size(x));
s= max(size(n));
si = max(size(x));

    for i = 1:s
        for z = 1:s
            Y(1,1) = Y(1,1) + (Interp(n, n(i,1), x(1,1))* Interp3(n, n(z,2), x(2,1))*fun(n(i,1),n(z,2)));
        end
    end

end




