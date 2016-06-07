function Y = funApprox(n, x )
%FUN~ Summary of this function goes here
%   Detailed explanation goes here
Y=0;
s= max(size(n));
t= max(size(x));
    for i=1:t
        for j = 1:s
            Y = Y +(Interp(n(:,i), n(j,i), x(1,i))*fun(n(j,i),0.5) );
        end
     
    end
end