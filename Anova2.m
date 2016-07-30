function Y = Anova2( n, x, w)
%UNTITLED Second Order Solution of Multi-Variate Problem. 
% this function is the 2 order approximation of  n variable function.

d=length(x);
Y=0;
value = 0;
value1 = 0;
   for k = 1:d
       E(k) = mean(n(:,k));
   end                                                             
for k= 1:d

a=unique(n(:,k));
     for l=(k+1):d
      b=unique(n(:,l));
      
      value = 0;
      value1 = 0;
      if abs(RelMean2(n, n(:,k),n(:,l), w, k, l, x, E))>10^(-5)
       c = length(a);
       t = length(b);
        for i = 1:c
            for j = 1:t
            e=E;
            e(k)=a(i);
            e(l)=b(j);
            Y = Y + (Interp(a, a(i), x(k))* Interp(b, b(j), x(l))*fun(e));
            end
        end
     
        for j = 1:t
        e=E;
        e(l)=b(j);
        v1 = Interp(b, b(j), x(l))*fun(e);
        value1 = value1 + v1;
        end
     
        for i = 1:c
        e=E;
        e(k)=a(i);
        v = Interp(a, a(i), x(k))*fun(e);
        value = value + v;
        end           
        Y=Y +fun(E)- value - value1;
     
     
      end
     end
     
       
end     
Y=Y+FApprox1(n,x,w);

