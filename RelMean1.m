function Value = RelMean1(n ,w, k, E )
% 
% 

value1 = 0;
value2 = 0;
value3 = 0;
a=unique(n);
c=length(n);
d=length(a);
for j=1:c 
value1 = 0;
    for i = 1:d
        e=E;
        e(k)=a(i);
        v = Interp( n, a(i), n(j))*fun(e);
        value1 = value1 + v;
    end
value2 = value2+((value1 - fun(E))*w(j)/c);
end
for i = 1:c
    value3 = value3 + fun(E)*w(i)/c;
end    
Value=value2/value3;
