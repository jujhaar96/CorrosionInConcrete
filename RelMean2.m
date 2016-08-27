function Value = RelMean2(n, N ,M ,w, k, K, x, E )
% 
% 

value1 = 0;
value2 = 0;
value3 = 0;
value4 = 0;
value5 = 0;
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

for i = 1:c
    X = [n(i,1), n(i,2)];
    value4 = value4 + FApprox1(n, X, w)*w(i)/c;
end  


Value= value3/value4;
end