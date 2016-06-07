function Y = FApprox2( n, x)
%UNTITLED Summary of this function goes here
%   j and x should be of the same dimension. we haven't provided the check
%   here. LATER.
Y=0;
a=unique(n(:,1));
b=unique(n(:,2));

s = length(unique(n));
    for i = 1:s
        for j = 1:s
            Y = Y + (Interp(a(:,1), a(i,1), x(1,1))* Interp(b(:,1), b(j,1), x(1,2))*fun(a(i,1),b(j,1)));
        end
    end

end




