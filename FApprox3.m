function Y = FApprox3( n, x,w)
%UNTITLED Second Order Solution of Multi-Variate Problem. 
% this function is the 2 order approximation of  n variable function.
Y=0;
d=length(x);
for k= 1:d
    a=unique(n(:,k));
    for l=(k+1):d
    b=unique(n(:,l));
        for m=(l+1):d;
            f=unique(n(:,m));
            c = length(a);
            t = length(b);
            z = length(f);
            for i = 1:c
                for j = 1:t
                    for u = 1:z 
                    e=ones(size(x))*0.5;
                    e(k)=a(i);
                    e(l)=b(j);
                    e(m)=f(u);
                    Y = Y + (Interp(a(:,1), a(i,1), x(k))*(Interp(f(:,1), f(u,1), x(m)))* Interp(b(:,1), b(j,1), x(l))*fun(e));
                    end
                end

            end
       end
    end
end
    f=(d-3)*FApprox2(n,x,w)
    q=((d-3)*(d-2)/2)*FApprox1(n,x,w)
    g=((d-1)*(d-2)*(d-3)/6)*fun(ones(size(x))*0.5)
    Y=Y-f-q-g;


