function [w, dwdx, dwdy,dwdz] = Weight33D(type, para, x,y,z,xI,yI,zI,dmI)
% EVALUATE WEIGHT FUNCTION
%
% SYNTAX: [w, dwdx, dwdy,dwdz] = Weight33D(type, para, x,y,z,xI,yI,zI,dmI)
%
% INPUT PARAMETERS
%    type - Type of weight function
%    para - Weight function parameter
%    x,y,z   - gauss point coordinates 
%    xI,yI,zI  -  nodal point coordinate
%    dmI - Support size
% OUTPUT PARAMETERS
%    w    - Value of weight function at r
%    dwdx - Value of first order derivative of weight function with respect to x at r
%    dwdy - Value of first order derivative of weight function with respect to y at r
%    dwdz - Value of first order derivative of weight function with respect to z at r
r1 = sqrt((x-xI).^2+(y-yI).^2+(z-zI).^2);
r  = r1/ dmI;   %   define the support size is a circle


        if(r==0)
            drdx=0;
            drdy=0;
            drdz=0;
        else
           drdx = x/(dmI.^2*r);
           drdy=y/(dmI.^2*r);
           drdz=z/(dmI.^2*r);
        end
        
        
% EVALUATE WEIGHT FUNCTION AND ITS FIRST AND SECOND ORDER OF DERIVATIVES WITH RESPECT r AT r

if (type == 'GAUSS')
   [w,dwdr] = Gauss(para,r);
elseif (type == 'CUBIC')
   [w,dwdr] = Cubic(r);
elseif (type == 'SPLI3')
   [w,dwdr] = Spline3(r);
elseif (type == 'SPLI5')
   [w,dwdr] = Spline5(r);
elseif (type == 'SPLIB')
   [w,dwdr] = BSpline(dmI/2,r);
elseif (type == 'power')
   [w,dwdr] = power_function(para,r);
elseif (type == 'CRBF1')
   [w,dwdr] = CSRBF1(r);
elseif (type == 'CRBF2')
   [w,dwdr] = CSRBF2(r);
elseif (type == 'CRBF3')
   [w,dwdr] = CSRBF3(r);
elseif (type == 'CRBF4')
   [w,dwdr] = CSRBF4(r);
elseif (type == 'CRBF5')
   [w,dwdr] = CSRBF5(r);
elseif (type == 'CRBF6')
   [w,dwdr] = CSRBF6(r);
elseif (type == 'LSQ11')
   [w,dwdr] = LSQ(r); %change to least square!
else
   error('Invalid type of weight function.');
end
   

dwdx  = dwdr * drdx;
dwdy= dwdr * drdy;
dwdz=dwdr*drdz;

function [w,dwdr] = Gauss(beta,r)
if (r>1.0)
   w     = 0.0;
   dwdr  = 0.0;
  else
   b2 = beta*beta;
   r2 = r*r;
   eb2 = exp(-b2);

   w     = (exp(-b2*r2) - eb2) / (1.0 - eb2);
   dwdr  = -2*b2*r*exp(-b2*r2) / (1.0 - eb2);
   
end

function [w,dwdr] = Cubic(r)
if (r>1.0)
   w     = 0.0;
   dwdr  = 0.0;
 
else
   w     = 1-6*r^2+8*r^3-3*r^4;
   dwdr  = -12*r+24*r^2-12*r^3;
  
end

function [w,dwdr] = Spline3(r)
if (r>1.0)
   w     = 0.0;
   dwdr  = 0.0;
 
elseif (r<=0.5)
   w     = 2/3 - 4*r^2 + 4*r^3;
   dwdr  = -8*r + 12*r^2;
 
else
   w     = 4/3 - 4*r + 4*r^2 - 4*r^3/3;
   dwdr  = -4 + 8*r -4*r^2;
  
end
function [w,dwdr] = Spline5(r)
if (r>1.0)
   w     = 0.0;
   dwdr  = 0.0;
 else
   w     = 1-10*r^3+15*r^4-6*r^5;
   dwdr  = -30*r^2 + 60*r^3-30*r^4;
  
end
function [w,dwdr] = BSpline(h,r)%h is the distance between two nodes.
if (r>1.0)
   w     = 0.0;
   dwdr  = 0.0;
 
elseif (r<=0.5)
   w     = 1/(pi*h^3)*(1-6*r^2+6*r^3);
   dwdr  = 1/(pi*h^3)*(-12*r+18*r^2);
 
else
   w     = 2/(pi*h^3)*(1-r)^3;
   dwdr  = -6/(pi*h^3)*(1-r)^2;
  
end


function [w,dwdr] = power_function(arfa,r)
if (r>1.0)
   w     = 0.0;
   dwdr  = 0.0;
else
    a2 = arfa*arfa;
   r2 = r*r;
    w     = exp(-r2/a2);
   dwdr  = (-2*r/a2)*exp(-r2/a2);
  
end

function [w,dwdr] = CSRBF2(r)
if (r>1.0)
   w     = 0.0;
   dwdr  = 0.0;
   
else
	w     = (1-r)^6*(6+36*r+82*r^2+72*r^3+30*r^4+5*r^5);
	dwdr  = 11*r*(r+2)*(5*r^3+15*r^2+18*r+4)*(r-1)^5;

end
function [w,dwdr] = CSRBF1(r)
if (r>1.0)
   w     = 0.0;
   dwdr  = 0.0;
   
else
	w     = (1-r)^4*(4+16*r+12*r^2+3*r^3);
	dwdr  = -4*(1-r)^3*(4+16*r+12*r^2+3*r^3)+(1-r)^4*(16+24*r+9*r^2);

end


function [w,dwdr] = CSRBF5(r)
if (r>1.0)
   w     = 0.0;
   dwdr  = 0.0;
   
else
	w     = (1-r)^6*(35*r^2+18*r+3);
	dwdr  =-6*(1-r)^5*(35*r^2+18*r+3)+(1-r)^6*(70*r+18);
   
end

function [w,dwdr] = CSRBF6(r)
if (r>1.0)
   w     = 0.0;
   dwdr  = 0.0;
   
else
	w     = (1-r)^8*(32*r^3+25*r^2+8*r+1);
	dwdr  =-8*(1-r)^7*(32*r^3+25*r^2+8*r+1)+(1-r)^8*(96*r^2+50*r+8);
   
end
function [w,dwdr] = LSQ(r)
if (r>1.0)
   w     = 0.0;
   dwdr  = 0.0;
   
else
	w     = 1;
	dwdr  =0;
end
function [w,dwdr] = CSRBF3(r)
if (r>1.0)
   w     = 0.0;
   dwdr  = 0.0;
elseif(r==0)
    w =  1/3;
    dwdr  =  0.0;
else
	w     = 1/3+r^2-4/3*r^3+2*r^2*log(r);
	dwdr  = 4*r-4*r^2+4*r*log(r);

end
function [w,dwdr] = CSRBF4(r)
if (r>1.0)
   w     = 0.0;
   dwdr  = 0.0; 
elseif(r==0)
    w =  1/15;
    dwdr  =  0.0;
   
else
	w     = 1/15+19/6*r^2-16/3*r^3+3*r^4-16/15*r^5+1/6*r^6+2*r^2*log(r);
	dwdr  = 25/3*r-16*r^2+12*r^3-16/3*r^4+r^5+4*r*log(r);

end