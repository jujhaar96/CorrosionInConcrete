function value = FApprox1(n, x, w )
%APPROXF First Order Solution of Multi-Variate Problem. 
% this function is the 1 order approximation of  n variable function.
value = 0;
d = length(x);
for j = 1:d
    E(j) = mean(n(:,j));
end
    for j = 1:d  
        a = unique(n(:,j));
        c = length(a);
        if abs(RelMean1(n(:,j),w,j, E ))>10^(-19)
            for i = 1:c
            e=E;
            e(j)=a(i);
            v = Interp(a, a(i), x(j))*fun(e);
            value = value+ v;
            end
            value = value - fun(E);      
        end
    end
value = value + fun(E);
end

