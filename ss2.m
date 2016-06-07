function Y = ss2( n, x)
%UNTITLED Summary of this function goes here
%   j and x should be of the same dimension. we haven't provided the check
%   here. LATER.
Y=0;%v
s= max(size(n));%n
    for i = 1:s % iterates through 1st column
    
         for z = 1:s % 2nd
           w = ones(size(x));
             for c = 1:s 
                for j = 1:s
                   if c~=i && j~=z
                   w = ((x(1,1)- n(c,1))./(n(i,1)-n(c,1))).*((x(1,2)- n(j,2))./(n(z,2)-n(j,2))).*w;
                   
                   end
                end
                 
             end
              Y = Y + w*fun(n(i,1),n(z,2))
         end
    
    end
end


