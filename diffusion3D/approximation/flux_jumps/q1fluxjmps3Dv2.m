function [jmp1,jmp2,jmp3,jmp4] = q1fluxjmps3Dv2(q1sol,efx,xyz,ev,ebound3D,s,t,l)
%Q1FLUXJMPS3DV2 vectorised flux jumps for element mid face node
%and edge nodes
%   jmp = q1fluxjmps3Dv2(q1sol,efx,xy,ev,ebound,s,t,l);
%   input
%          q1sol        vertex solution vector
%          efx          presorted element connectivity array
%          xyz          vertex coordinate vector
%          ev           element mapping matrix
%          ebound3d     element edge boundary matrix
%          s            Gaussian point , s\in (-1,1)
%          t            Gaussian point , t\in (-1,1)
%          l            Gaussian point , l\in (-1,1)
%
%   output
%          jmp1          component elementwise edge flux jumps
%          jmp2          component elementwise edge flux jumps
%          jmp3          component elementwise edge flux jumps
%          jmp4          component elementwise edge flux jumps
%
% IFISS function: GP 22 June 2022; .
% Copyright (c)   G.Papanikos, C.E. Powell, D.J. Silvester
tic
x=xyz(:,1); y=xyz(:,2); z=xyz(:,3); nvtx=length(x);
nel=length(ev(:,1));

% initialise global matrices
jmp1   = zeros(nel,6);
flux1  = zeros(nel,6);
jmp2   = zeros(nel,6);
flux2  = zeros(nel,6);
jmp3   = zeros(nel,6);
flux3  = zeros(nel,6);
jmp4   = zeros(nel,6);
flux4  = zeros(nel,6);

zero_v=zeros(nel,1); one_v=ones(nel,1);
%
% inner loop over elements
xl_v = x(ev);
yl_v = y(ev);
zl_v = z(ev);
sl_v = q1sol(ev);
% evaluate derivatives and normal flux on a face in turn
% and assemble into flux and face area adjacency matrices

%% face 1 nodes

% ns= s, nt=-1 , nl =-1
[~,invjac,~,~,~,dphidz]  = deriv3D(s,-1,-1,xl_v,yl_v,zl_v);
fx_v=zero_v;
for  ivtx=1:8
    fx_v = fx_v +(-dphidz(:,ivtx)).*invjac(:).*sl_v(:,ivtx);
end
flux_node9_face1=fx_v;

% ns= 1, nt=t , nl =-1
[~,invjac,~,~,~,dphidz]  = deriv3D(1,t,-1,xl_v,yl_v,zl_v);
fx_v=zero_v;
for  ivtx=1:8
    fx_v = fx_v +(-dphidz(:,ivtx)).*invjac(:).*sl_v(:,ivtx);
end
flux_node15_face1=fx_v;

% ns= 0, nt=1 , nl =-1
[~,invjac,~,~,~,dphidz]  = deriv3D(s,1,-1,xl_v,yl_v,zl_v);
fx_v=zero_v;
for  ivtx=1:8
    fx_v = fx_v +(-dphidz(:,ivtx)).*invjac(:).*sl_v(:,ivtx);
end
flux_node23_face1=fx_v;

% ns= -1, nt=0 , nl =-1
[~,invjac,~,~,~,dphidz]  = deriv3D(-1,t,-1,xl_v,yl_v,zl_v);
fx_v=zero_v;
for  ivtx=1:8
    fx_v = fx_v +(-dphidz(:,ivtx)).*invjac(:).*sl_v(:,ivtx);
end
flux_node21_face1=fx_v;

flux1(:,1) = flux_node9_face1;
flux2(:,1) = flux_node15_face1;
flux3(:,1) = flux_node23_face1;
flux4(:,1) = flux_node21_face1;
%% face 2 nodes

% ns = 1, nt= 0 , nl= -1
[~,invjac,~,dphidx,~,~]  = deriv3D(1,t,-1,xl_v,yl_v,zl_v);
fx_v=zero_v;
for  ivtx=1:8
    fx_v = fx_v +(dphidx(:,ivtx)).*invjac(:).*sl_v(:,ivtx);
end
flux_node15_face2=fx_v;

% ns = 1, nt= 1 , nl= 0
[~,invjac,~,dphidx,~,~]  = deriv3D(1,1,l,xl_v,yl_v,zl_v);
fx_v=zero_v;
for  ivtx=1:8
    fx_v = fx_v +(dphidx(:,ivtx)).*invjac(:).*sl_v(:,ivtx);
end
flux_node24_face2=fx_v;

% ns = 1, nt= 0 , nl= 1
[~,invjac,~,dphidx,~,~]  = deriv3D(1,t,1,xl_v,yl_v,zl_v);
fx_v=zero_v;
for  ivtx=1:8
    fx_v = fx_v +(dphidx(:,ivtx)).*invjac(:).*sl_v(:,ivtx);
end
flux_node17_face2=fx_v;

% ns = 1, nt= -1 , nl= 0
[~,invjac,~,dphidx,~,~]  = deriv3D(1,-1,l,xl_v,yl_v,zl_v);
fx_v=zero_v;
for  ivtx=1:8
    fx_v = fx_v +(dphidx(:,ivtx)).*invjac(:).*sl_v(:,ivtx);
end
flux_node10_face2=fx_v;

flux1(:,2) = flux_node15_face2;
flux2(:,2) = flux_node24_face2;
flux3(:,2) = flux_node17_face2;
flux4(:,2) = flux_node10_face2;

%% face 3 nodes
% ns= 0, nt= -1 nl=1
[~,invjac,~,~,~,dphidz]  = deriv3D(s,-1,1,xl_v,yl_v,zl_v);
fx_v=zero_v;
for  ivtx=1:8
    fx_v = fx_v +( dphidz(:,ivtx)).*invjac(:).*sl_v(:,ivtx);
end
flux_node11_face3=fx_v;

% ns= 1, nt= 0 nl=1
[~,invjac,~,~,~,dphidz]  = deriv3D(1,t,1,xl_v,yl_v,zl_v);
fx_v=zero_v;
for  ivtx=1:8
    fx_v = fx_v +( dphidz(:,ivtx)).*invjac(:).*sl_v(:,ivtx);
end
flux_node17_face3=fx_v;

% ns= 0, nt= 1 nl=1
[~,invjac,~,~,~,dphidz]  = deriv3D(s,1,1,xl_v,yl_v,zl_v);
fx_v=zero_v;
for  ivtx=1:8
    fx_v = fx_v +( dphidz(:,ivtx)).*invjac(:).*sl_v(:,ivtx);
end
flux_node25_face3=fx_v;

% ns= -1, nt= 0 nl=1
[~,invjac,~,~,~,dphidz]  = deriv3D(-1,t,1,xl_v,yl_v,zl_v);
fx_v=zero_v;
for  ivtx=1:8
    fx_v = fx_v +( dphidz(:,ivtx)).*invjac(:).*sl_v(:,ivtx);
end
flux_node19_face3=fx_v;

flux1(:,3) = flux_node11_face3;
flux2(:,3) = flux_node17_face3;
flux3(:,3) = flux_node25_face3;
flux4(:,3) = flux_node19_face3;

%% face 4 nodes
% ns = -1, nt= 0, nl =-1
[~,invjac,~,dphidx,~,~]  = deriv3D(-1,t,-1,xl_v,yl_v,zl_v);
fx_v=zero_v;
for  ivtx=1:8
    fx_v = fx_v +(-dphidx(:,ivtx)).*invjac(:).*sl_v(:,ivtx);
end
flux_node21_face4=fx_v;
% ns = -1, nt= 1, nl =0
[~,invjac,~,dphidx,~,~]  = deriv3D(-1,1,l,xl_v,yl_v,zl_v);
fx_v=zero_v;
for  ivtx=1:8
    fx_v = fx_v +(-dphidx(:,ivtx)).*invjac(:).*sl_v(:,ivtx);
end
flux_node26_face4=fx_v;
% ns = -1, nt= 0, nl =1
[~,invjac,~,dphidx,~,~]  = deriv3D(-1,t,1,xl_v,yl_v,zl_v);
fx_v=zero_v;
for  ivtx=1:8
    fx_v = fx_v +(-dphidx(:,ivtx)).*invjac(:).*sl_v(:,ivtx);
end
flux_node19_face4=fx_v;
% ns = -1, nt= -1, nl =0
[~,invjac,~,dphidx,~,~]  = deriv3D(-1,-1,l,xl_v,yl_v,zl_v);
fx_v=zero_v;
for  ivtx=1:8
    fx_v = fx_v +(-dphidx(:,ivtx)).*invjac(:).*sl_v(:,ivtx);
end
flux_node12_face4=fx_v;

flux1(:,4) = flux_node21_face4;
flux2(:,4) = flux_node26_face4;
flux3(:,4) = flux_node19_face4;
flux4(:,4) = flux_node12_face4;

%% face 5 nodes
% ns=0, nt=-1 , nl= -1
[~,invjac,~,~,dphidy,~]  = deriv3D(s,-1,-1,xl_v,yl_v,zl_v);
fx_v=zero_v;
for  ivtx=1:8
    fx_v = fx_v +(-dphidy(:,ivtx)).*invjac(:).*sl_v(:,ivtx);
end
flux_node9_face5=fx_v;

% ns=1, nt=-1 , nl= 0
[~,invjac,~,~,dphidy,~]  = deriv3D(1,-1,l,xl_v,yl_v,zl_v);
fx_v=zero_v;
for  ivtx=1:8
    fx_v = fx_v +(-dphidy(:,ivtx)).*invjac(:).*sl_v(:,ivtx);
end
flux_node10_face5=fx_v;

% ns=0, nt=-1 , nl= 1
[~,invjac,~,~,dphidy,~]  = deriv3D(s,-1,1,xl_v,yl_v,zl_v);
fx_v=zero_v;
for  ivtx=1:8
    fx_v = fx_v +(-dphidy(:,ivtx)).*invjac(:).*sl_v(:,ivtx);
end
flux_node11_face5=fx_v;

% ns=-1, nt=-1 , nl= 0
[~,invjac,~,~,dphidy,~]  = deriv3D(-1,-1,l,xl_v,yl_v,zl_v);
fx_v=zero_v;
for  ivtx=1:8
    fx_v = fx_v +(-dphidy(:,ivtx)).*invjac(:).*sl_v(:,ivtx);
end
flux_node12_face5=fx_v;

flux1(:,5) = flux_node9_face5;
flux2(:,5) = flux_node10_face5;
flux3(:,5) = flux_node11_face5;
flux4(:,5) = flux_node12_face5;


%% face 6 nodes
% ns= 0, nt= 1, nl =-1
[~,invjac,~,~,dphidy,~]  = deriv3D(s,1,-1,xl_v,yl_v,zl_v);
fx_v=zero_v;
for  ivtx=1:8
    fx_v = fx_v +(dphidy(:,ivtx)).*invjac(:).*sl_v(:,ivtx);
end
flux_node23_face6=fx_v;

% ns= 1, nt= 1, nl =0
[~,invjac,~,~,dphidy,~]  = deriv3D(1,1,l,xl_v,yl_v,zl_v);
fx_v=zero_v;
for  ivtx=1:8
    fx_v = fx_v +(dphidy(:,ivtx)).*invjac(:).*sl_v(:,ivtx);
end
flux_node24_face6=fx_v;

% ns= 0, nt= 1, nl =1
[~,invjac,~,~,dphidy,~]  = deriv3D(s,1,1,xl_v,yl_v,zl_v);
fx_v=zero_v;
for  ivtx=1:8
    fx_v = fx_v +(dphidy(:,ivtx)).*invjac(:).*sl_v(:,ivtx);
end
flux_node25_face6=fx_v;

% ns= -1, nt= 1, nl =0
[~,invjac,~,~,dphidy,~]  = deriv3D(-1,1,l,xl_v,yl_v,zl_v);
fx_v=zero_v;
for  ivtx=1:8
    fx_v = fx_v +(dphidy(:,ivtx)).*invjac(:).*sl_v(:,ivtx);
end
flux_node26_face6=fx_v;

flux1(:,6) = flux_node23_face6;
flux2(:,6) = flux_node24_face6;
flux3(:,6) = flux_node25_face6;
flux4(:,6) = flux_node26_face6;

%
%% evaluate flux jump on each face in turn

jmp1(1:nel,1) = flux1(1:nel,1)+flux1(efx(1:nel,1),3);
jmp2(1:nel,1) = flux2(1:nel,1)+flux2(efx(1:nel,1),3);
jmp3(1:nel,1) = flux3(1:nel,1)+flux3(efx(1:nel,1),3);
jmp4(1:nel,1) = flux4(1:nel,1)+flux4(efx(1:nel,1),3);
%
jmp1(1:nel,2) = flux1(1:nel,2)+flux1(efx(1:nel,2),4);
jmp2(1:nel,2) = flux2(1:nel,2)+flux2(efx(1:nel,2),4);
jmp3(1:nel,2) = flux3(1:nel,2)+flux3(efx(1:nel,2),4);
jmp4(1:nel,2) = flux4(1:nel,2)+flux4(efx(1:nel,2),4);
%
jmp1(1:nel,3) = flux1(1:nel,3)+flux1(efx(1:nel,3),1);
jmp2(1:nel,3) = flux2(1:nel,3)+flux2(efx(1:nel,3),1);
jmp3(1:nel,3) = flux3(1:nel,3)+flux3(efx(1:nel,3),1);
jmp4(1:nel,3) = flux4(1:nel,3)+flux4(efx(1:nel,3),1);
%
jmp1(1:nel,4) = flux1(1:nel,4)+flux1(efx(1:nel,4),2);
jmp2(1:nel,4) = flux2(1:nel,4)+flux2(efx(1:nel,4),2);
jmp3(1:nel,4) = flux3(1:nel,4)+flux3(efx(1:nel,4),2);
jmp4(1:nel,4) = flux4(1:nel,4)+flux4(efx(1:nel,4),2);
%
jmp1(1:nel,5) = flux1(1:nel,5)+flux1(efx(1:nel,5),6);
jmp2(1:nel,5) = flux2(1:nel,5)+flux2(efx(1:nel,5),6);
jmp3(1:nel,5) = flux3(1:nel,5)+flux3(efx(1:nel,5),6);
jmp4(1:nel,5) = flux4(1:nel,5)+flux4(efx(1:nel,5),6);
%
jmp1(1:nel,6) = flux1(1:nel,6)+flux1(efx(1:nel,6),5);
jmp2(1:nel,6) = flux2(1:nel,6)+flux2(efx(1:nel,6),5);
jmp3(1:nel,6) = flux3(1:nel,6)+flux3(efx(1:nel,6),5);
jmp4(1:nel,6) = flux4(1:nel,6)+flux4(efx(1:nel,6),5);
%
etime=toc;
%  fprintf('flux jump assembly took %6.3e seconds\n',etime)
%
% remove Dirichlet boundary edge contributions
nbde=length(ebound3D(:,1));
for k=1:nbde
    jmp1(ebound3D(k,1),ebound3D(k,2))=0;
    jmp2(ebound3D(k,1),ebound3D(k,2))=0;
    jmp3(ebound3D(k,1),ebound3D(k,2))=0;
    jmp4(ebound3D(k,1),ebound3D(k,2))=0;
end
%
return

