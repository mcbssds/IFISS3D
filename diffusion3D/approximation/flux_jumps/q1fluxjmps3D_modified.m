function jmp = q1fluxjmps3D_modified(q1sol,efx,xyz,ev,ebound3D,s,t)
%Q1FLUXJMPS3D_MODIFIED vectorised flux jumps on the midpoint and edges of
%the faces for hexahedral element
%   jmp = q1fluxjmps3D_modified(q1sol,efx,xyz,ev,ebound3D,s,t)
%   input
%          q1sol        vertex solution vector
%          efx          presorted element connectivity array
%          xyz          vertex coordinate vector
%          ev           element mapping matrix
%          ebound3d     element edge boundary matrix
%          s            quadrature point coordinate in (-1,1)
%          t            quadrature point coordinate in (-1,1)
%   output
%          jmp          component elementwise edge flux jumps
%
% call function  deriv3D 
% IFISS function: GP 9 June 2022 
% Copyright (c) 2022  G.Papanikos, C.E. Powell, D.J. Silvester
tic
x=xyz(:,1); y=xyz(:,2); z=xyz(:,3); nvtx=length(x);
nel=length(ev(:,1));

% initialise global matrices
jmp   = zeros(nel,6);
flux  = zeros(nel,6);
zero_v=zeros(nel,1); one_v=ones(nel,1);
%

xl_v = x(ev);
yl_v = y(ev);
zl_v = z(ev);
sl_v = q1sol(ev);
% evaluate derivatives and normal flux on a face in turn
% and assemble into flux and face area adjacency matrices

% ns= 0, nt=0 , nl =-1
[~,invjac,~,~,~,dphidz]  = deriv3D(s,t,-1,xl_v,yl_v,zl_v);
fx_v = sum((-dphidz).*invjac(:).*sl_v,2);
flux(:,1)=fx_v;
%
% ns = 1, nt= 0 , nl= 0
[~,invjac,~,dphidx,~,~]  = deriv3D(1,s,t,xl_v,yl_v,zl_v);
fx_v = sum(dphidx.*invjac(:).*sl_v,2);
flux(:,2)=fx_v;
%
% ns= 0, nt= 0 nl=1
[~,invjac,~,~,~,dphidz]  = deriv3D(s,t,1,xl_v,yl_v,zl_v);
fx_v = sum(dphidz.*invjac(:).*sl_v,2);
flux(:,3)=fx_v;
%
% ns = -1, nt= 0, nl =0
[~,invjac,~,dphidx,~,~]  = deriv3D(-1,s,t,xl_v,yl_v,zl_v);
fx_v = sum(-dphidx.*invjac(:).*sl_v,2);
flux(:,4)=fx_v;

% ns=0, nt=-1 , nl= 0
[~,invjac,~,~,dphidy,~]  = deriv3D(s,-1,t,xl_v,yl_v,zl_v);
fx_v = sum(-dphidy.*invjac(:).*sl_v,2);
flux(:,5)=fx_v;
%

% ns= 0, nt= 1, nl =0
[~,invjac,~,~,dphidy,~]  = deriv3D(s,1,t,xl_v,yl_v,zl_v);
fx_v = sum(dphidy.*invjac(:).*sl_v,2);
flux(:,6)=fx_v;

%
% evaluate flux jump on each face in turn
% nz= -1, ny= 0, nz =0
jmp(1:nel,1) = flux(1:nel,1)+flux(efx(1:nel,1),3);
% nx= 1, ny=0, nz=0
jmp(1:nel,2) = flux(1:nel,2)+flux(efx(1:nel,2),4);
%
% nx= 0, ny=0, nz=1
jmp(1:nel,3) = flux(1:nel,3)+flux(efx(1:nel,3),1);
%
% nx= -1, ny=0 nz =0
jmp(1:nel,4) = flux(1:nel,4)+flux(efx(1:nel,4),2);
%
% nx= 0, ny=-1 nz = 0
jmp(1:nel,5) = flux(1:nel,5)+flux(efx(1:nel,5),6);
%
% nx= 0, ny=1, nz = 0
jmp(1:nel,6) = flux(1:nel,6)+flux(efx(1:nel,6),5);
%
etime=toc;
%  fprintf('flux jump assembly took %6.3e seconds\n',etime)
%
% remove Dirichlet boundary face contributions
nbde=length(ebound3D(:,1));
for k=1:nbde
    jmp(ebound3D(k,1),ebound3D(k,2))=0;
end
%
return
