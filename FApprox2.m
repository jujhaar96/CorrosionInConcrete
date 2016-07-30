function Y = FApprox2( n, x,w)
%UNTITLED Second Order Solution of Multi-Variate Problem. 
% this function is the 2 order approximation of  n variable function.
Y=0;
d=length(x);
for k= 1:d
    a=unique(n(:,k));
    for l=(k+1):d
    b=unique(n(:,l));
    c = length(a);
    t = length(b);
        for i = 1:c
            for j = 1:t
            e=ones(size(x))*0.5;
            e(k)=a(i);
            e(l)=b(j);
            Y = Y + (Interp(a(:,1), a(i,1), x(k))* Interp(b(:,1), b(j,1), x(l))*fun(e));
            end
        end
    end
end

    f=(d-2)*FApprox1(n,x,w);
    q=(d-1)*(d-2)*0.5*fun(ones(size(x))*0.5);
    Y=Y-f-q;


