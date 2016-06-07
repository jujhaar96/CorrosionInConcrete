function Y = Elements(n)

a = length(n);
Y = n;

index=0;
for i = 1:(a-1)
     
     for j = (i+1):a
         
         if (max(abs((n(i)-n(j)))) == 0)
             index = [index,j];
         end
     end
end
index
s = length(index);
i=s;
while i>1
    
Y(index(i))=[];
i=i-1;
end
