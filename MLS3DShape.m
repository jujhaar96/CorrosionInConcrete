function [PHI, DPHIx, DPHIy,DPHIz] = MLS3DShape(m, nnodes, xI,yI,zI, npoints, x,y, z,dmI, type, para)
% SHAPE FUNCTION OF 3D MLS APPROXIMATION
%
% SYNTAX: [PHI, DPHIx, DPHIy,DPHIz] = MLS3DShape(m, nnodes, xI,yI,zI, npoints, x,y, z,dmI, type, para)
%
% INPUT PARAMETERS
%    m - Total number of basis functions (1: Constant basis;  2: Linear basis;  3: Quadratic basis)
%    nnodes  - Total number of nodes used to construct MLS approximation
%    npoints - Total number of points whose MLS shape function to be evaluated
%    xI,yI,zI(nnodes) - Coordinates of nodes used to construct MLS approximation
%    xi,yi,zi(npoints) - Coordinates of points whose MLS shape function to be evaluated
%    dm(nnodes) - Radius of support of nodes
%    type - Type of weight function
%    para  - Weight function parameter
%
% OUTPUT PARAMETERS
%    PHI   - MLS Shpae function
%    DPHIx - First order derivatives of MLS Shpae function to x
%    DPHIy - First order derivatives of MLS Shpae function to y
%    DPHIz - First order derivatives of MLS Shpae function to z
% INITIALIZE WEIGHT FUNCTION MATRICES
DmI = [];
wI   = zeros (1, nnodes);  % Weight funciton
dwdxI  = zeros (1, nnodes);
dwdyI = zeros (1, nnodes);
dwdzI = zeros (1, nnodes);
xII = zeros(1,nnodes);
yII = zeros(1,nnodes);
zII = zeros(1,nnodes);
% INITIALIZE SHAPE FUNCTION MATRICES
PHI   = zeros(npoints, nnodes);
DPHIx  = zeros(npoints, nnodes);
DPHIy = zeros(npoints, nnodes);
DPHIz = zeros(npoints, nnodes);
% LOOP OVER ALL EVALUATION POINTS TO CALCULATE VALUE OF SHAPE FUNCTION Fi(X)
for j = 1 : npoints
    DmI = dmI;
	% DETERMINE WEIGHT FUNCTIONS AND THEIR DERIVATIVES AT EVERY NODE
	for i = 1 : nnodes
		[wI(i), dwdxI(i), dwdyI(i), dwdzI(i)] =rectangleWeight3D(type, para, x(j),y(j),z(j),xI(i),yI(i),zI(i),DmI(i),DmI(i),DmI(i));
       xII(1,i)=xI(i);
       yII(1,i)=yI(i);
       zII(1,i)=zI(i);
	end
   
   % EVALUATE BASIS p, B MATRIX AND THEIR DERIVATIVES
   if (m == 1)  % Shepard function
      p = [ones(1, nnodes)]; 
      pxyz   = [1];
      dpdx  = [0];
      dpdy = [0];
      dpdz = [0];
      
      B    = p .* [wI];
      DBdx   = p .* [dwdxI];
      DBdy  = p .* [dwdyI];
      DBdz  = p .* [dwdzI];
   elseif (m == 4)
      p = [ones(1, nnodes); xII;yII;zII]; 
      pxyz   = [1; x(j);y(j);z(j)];
      dpdx  = [0; 1;0;0];
      dpdy = [0; 0;1;0];
      dpdz = [0; 0;0;1];
      
      B    = p .* [wI; wI;wI;wI];
      DBdx   = p .* [dwdxI; dwdxI;dwdxI;dwdxI];
      DBdy  = p .* [dwdyI; dwdyI;dwdyI;dwdyI];
      DBdz  = p .* [dwdzI; dwdzI;dwdzI;dwdzI];
   elseif (m == 10)
      p = [ones(1, nnodes); xII;yII; zII;xII.*xII;xII.*yII;yII.*yII;yII.*zII;zII.*zII;xII.*zII]; 
      pxyz   = [1; x(j); y(j);z(j);x(j)*x(j);x(j)*y(j);y(j)*y(j);y(j)*z(j);z(j)*z(j);x(j)*z(j)];
      dpdx  = [0; 1; 0;0;2*x(j);y(j);0;0;0;z(j)];
      dpdy = [0;0;1;0;0;x(j);2*y(j);z(j);0;0];
      dpdz =[0;0;0;1;0;0;0;y(j);2*z(j);x(j)];
      B    = p .* [wI; wI; wI; wI; wI; wI; wI; wI; wI;wI];
      DBdx   = p .* [dwdxI; dwdxI;dwdxI;dwdxI; dwdxI;dwdxI;dwdxI;dwdxI; dwdxI;dwdxI];
      DBdy  = p .* [dwdyI; dwdyI;dwdyI;dwdyI; dwdyI;dwdyI;dwdyI;dwdyI; dwdyI;dwdyI];
      DBdz  = p .* [dwdzI; dwdzI;dwdzI;dwdzI; dwdzI;dwdzI;dwdzI;dwdzI; dwdzI;dwdzI];
   else
      error('Invalid order of basis.');
   end
   
   % EVALUATE MATRICES A AND ITS DERIVATIVES
	A   = zeros (m, m);
	DAdx  = zeros (m, m);
	DAdy= zeros (m, m);
    DAdz= zeros (m, m);
for i = 1 : nnodes
      pp = p(:,i) * p(:,i)';
      
      A   = A   + wI(i) * pp;
      DAdx  = DAdx  + dwdxI(i) * pp;
      DAdy = DAdy + dwdyI(i) * pp;
      DAdz = DAdz + dwdzI(i) * pp;
end
   ARcond=  rcond(A);
   
while ARcond<=9.9999999e-015   % 判断条件数
    DmI=1.1*DmI;
    for i = 1 : nnodes
		[wI(i), dwdxI(i), dwdyI(i), dwdzI(i)] =rectangleWeight3D(type, para, x(j),y(j),z(j),xI(i),yI(i),zI(i),DmI(i),DmI(i),DmI(i));
       xII(1,i)=xI(i);
       yII(1,i)=yI(i);
       zII(1,i)=zI(i);
	end
   
   % EVALUATE BASIS p, B MATRIX AND THEIR DERIVATIVES
   if (m == 1)  % Shepard function
      p = [ones(1, nnodes)]; 
      pxyz   = [1];
      dpdx  = [0];
      dpdy = [0];
      dpdz = [0];
      
      B    = p .* [wI];
      DBdx   = p .* [dwdxI];
      DBdy  = p .* [dwdyI];
      DBdz  = p .* [dwdzI];
   elseif (m == 4)
      p = [ones(1, nnodes); xII;yII;zII]; 
      pxyz   = [1; x(j);y(j);z(j)];
      dpdx  = [0; 1;0;0];
      dpdy = [0; 0;1;0];
      dpdz = [0; 0;0;1];
      
      B    = p .* [wI; wI;wI;wI];
      DBdx   = p .* [dwdxI; dwdxI;dwdxI;dwdxI];
      DBdy  = p .* [dwdyI; dwdyI;dwdyI;dwdyI];
      DBdz  = p .* [dwdzI; dwdzI;dwdzI;dwdzI];
   elseif (m == 10)
      p = [ones(1, nnodes); xII;yII; zII;xII.*xII;xII.*yII;yII.*yII;yII.*zII;zII.*zII;xII.*zII]; 
      pxyz   = [1; x(j); y(j);z(j);x(j)*x(j);x(j)*y(j);y(j)*y(j);y(j)*z(j);z(j)*z(j);x(j)*z(j)];
      dpdx  = [0; 1; 0;0;2*x(j);y(j);0;0;0;z(j)];
      dpdy = [0;0;1;0;0;x(j);2*y(j);z(j);0;0];
      dpdz =[0;0;0;1;0;0;0;y(j);2*z(j);x(j)];
      B    = p .* [wI; wI; wI; wI; wI; wI; wI; wI; wI;wI];
      DBdx   = p .* [dwdxI; dwdxI;dwdxI;dwdxI; dwdxI;dwdxI;dwdxI;dwdxI; dwdxI;dwdxI];
      DBdy  = p .* [dwdyI; dwdyI;dwdyI;dwdyI; dwdyI;dwdyI;dwdyI;dwdyI; dwdyI;dwdyI];
      DBdz  = p .* [dwdzI; dwdzI;dwdzI;dwdzI; dwdzI;dwdzI;dwdzI;dwdzI; dwdzI;dwdzI];
   else
      error('Invalid order of basis.');
   end
   
   % EVALUATE MATRICES A AND ITS DERIVATIVES
	A   = zeros (m, m);
	DAdx  = zeros (m, m);
	DAdy= zeros (m, m);
    DAdz= zeros (m, m);
for i = 1 : nnodes
      pp = p(:,i) * p(:,i)';
      
      A   = A   + wI(i) * pp;
      DAdx  = DAdx  + dwdxI(i) * pp;
      DAdy = DAdy + dwdyI(i) * pp;
      DAdz = DAdz + dwdzI(i) * pp;
end
   ARcond=  rcond(A); 
   
end                %判断A矩阵的条件数


       AInv = inv(A);
      
   rxyz  = AInv * pxyz;
   PHI(j,:) = rxyz' * B;   % shape function
    
   drdx  = AInv * (dpdx -DAdx* rxyz);
   DPHIx(j,:) = drdx' * B + rxyz' * DBdx;   % first order derivatives of shape function with respect to x
   
     drdy = AInv * (dpdy -DAdy* rxyz);
   DPHIy(j,:) = drdy' * B + rxyz' * DBdy;  % first order derivatives of shape function to y
   
   drdz = AInv * (dpdz -DAdz* rxyz);    % first order derivatives of shape function to y
   DPHIz(j,:) = drdz' * B + rxyz' * DBdz;
  
end