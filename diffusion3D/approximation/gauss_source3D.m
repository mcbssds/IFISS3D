function ff = gauss_source3D(s,t,l,xl,yl,zl,ngpt)
%GAUSS_SOURCE3D evaluates source term at Gauss point
%   ff = gauss_source3D(s,t,xl,yl,zl);
%   inputs:
%          s         reference element x coordinate   
%          t         reference element y coordinate
%          l         reference element y coordinate
%          xl        physical element x vertex coordinates 
%          yl        physical element y vertex coordinates  
%          zl        physical element z vertex coordinates
%  output:
%          ff        vector of f values for all elements
% IFISS function: GP; 9 June 2022.
% Copyright(c) 2022  G.Papanikos,  C.E. Powell, D.J. Silvester

if nargin<7
   ngpt = 1;
end    
[nel,nv] = size(xl);
zero_v = zeros(nel,ngpt); xx=zero_v; yy=xx; zz=xx;
[phi_e,~,~,~] = shape3D(s,t,l);
 
for ivtx=1:8
    xx = xx + phi_e(ivtx) * xl(:,ivtx);
    yy = yy + phi_e(ivtx) * yl(:,ivtx);
    zz = zz + phi_e(ivtx) * zl(:,ivtx);
end

ff=specific_rhs3D(xx,yy,zz,nel,ngpt);
return
