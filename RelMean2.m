function Value = RelMean2(n, N ,M ,w, k, K, x, E )
% 
% 

value1 = 0;
value2 = 0;
value3 = 0;
a=unique(N);
b=unique(M);
c=length(N);
d=length(a);
f=length(b);
Y=0;

for y=1:c
  
 for z=(y+1):c
 
     Y=0;
     value1 = 0;
     value2 = 0;
     for i = 1:d
            for j = 1:f
            e=E;
            e(k)=a(i);
            e(K)=b(j);
            Y = Y + (Interp(a, a(i), N(y))* Interp(b, b(j), M(y))*fun(e));
            end
    end
     for i = 1:d
        e=E;
        e(k)=a(i);
        v = Interp( N, a(i), N(j))*fun(e);
        value1 = value1 + v;
     end
     for i = 1:f
        e=E;
        e(K)=a(i);
        v = Interp( N, a(i), N(j))*fun(e);
        value2 = value2 + v;
    end
    
    
    
    
  value3 = value3+((Y - value1 - value2 + fun(E))*w(y)/c);  
  end
end
Value= value3/(FApprox1(n,x,w));
end