function jmp = q1fluxjmps3D(q1sol,efx,xyz,ev,ebound3D,s)
%Q1FLUXJMPS3D vectorised flux jumps at midpoints of element faces
%   jmp = q1fluxjmps3D(q1sol,efx,xy,ev,ebound3D,s);
%   inputs:
%          q1sol        vertex solution vector
%          efx          presorted element connectivity array
%          xyz          vertex coordinate vector
%          ev           element mapping matrix
%          ebound3d     element edge boundary matrix
%          s            quadrature point coordinate in (-1,1)
%   outputs:
%          jmp          component elementwise edge flux jumps
%
% IFISS function: GP 23 June 2022 (modified from DJS; 27 September 2013).
% Copyright (c)   G.Papanikos, C.E. Powell, D.J. Silvester
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
[~,invjac,~,~,~,dphidz]  = deriv3D(s,s,-1,xl_v,yl_v,zl_v);
fx_v=zero_v;
for  ivtx=1:8
    fx_v = fx_v +(-dphidz(:,ivtx)).*invjac(:).*sl_v(:,ivtx);
end
flux(:,1)=fx_v;
%
% ns = 1, nt= 0 , nl= 0
[~,invjac,~,dphidx,~,~]  = deriv3D(1,s,s,xl_v,yl_v,zl_v);
fx_v=zero_v;
for  ivtx=1:8
    fx_v = fx_v +(dphidx(:,ivtx)).*invjac(:).*sl_v(:,ivtx);
end
flux(:,2)=fx_v;
%
% ns= 0, nt= 0 nl=1
[~,invjac,~,~,~,dphidz]  = deriv3D(s,s,1,xl_v,yl_v,zl_v);
fx_v=zero_v;
for  ivtx=1:8
    fx_v = fx_v +( dphidz(:,ivtx)).*invjac(:).*sl_v(:,ivtx);
end
flux(:,3)=fx_v;
%
% ns = -1, nt= 0, nl =0
[~,invjac,~,dphidx,~,~]  = deriv3D(-1,s,s,xl_v,yl_v,zl_v);
fx_v=zero_v;
for  ivtx=1:8
    fx_v = fx_v +(-dphidx(:,ivtx)).*invjac(:).*sl_v(:,ivtx);
end
flux(:,4)=fx_v;

% ns=0, nt=-1 , nl= 0
[~,invjac,~,~,dphidy,~]  = deriv3D(s,-1,s,xl_v,yl_v,zl_v);
fx_v=zero_v;
for  ivtx=1:8
    fx_v = fx_v +(-dphidy(:,ivtx)).*invjac(:).*sl_v(:,ivtx);
end
flux(:,5)=fx_v;
%

% ns= 0, nt= 1, nl =0
[~,invjac,~,~,dphidy,~]  = deriv3D(s,1,s,xl_v,yl_v,zl_v);
fx_v=zero_v;
for  ivtx=1:8
    fx_v = fx_v +(dphidy(:,ivtx)).*invjac(:).*sl_v(:,ivtx);
end
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
