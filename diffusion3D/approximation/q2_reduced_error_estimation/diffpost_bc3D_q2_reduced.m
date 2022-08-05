function [elerr_p] = diffpost_bc3D_q2_reduced(elerror,fez,xyz,ev,ebound3D);
%DIFFPOST_BC3D_Q2_REDUCED postprocesses local Poisson error estimator
%   [elerr_p] = diffpost_bc3D_q2_reduced(elerror,fez,xyz,ev,ebound3D);
%   inputs:
%          elerror       element error estimate (without BC imposed)
%          fez           elementwise rhs vectors
%          xyz           vertex coordinate vector
%          ev            element mapping matrix
%          ebound3D      element face boundary matrix
%   output:
%          elerr_p       element error estimate with BC correction
%
%   calls functions gausspoints_oned, gausspoints_threed deriv3D qderiv3D
%   IFISS function: GP; 09 June 2022
% Copyright (c)  2022 G.Papanikos, C.E. Powell, D.J. Silvester
x=xyz(:,1); y=xyz(:,2); z= xyz(:,3);
nel=length(ev(:,1));
lev=[ev,ev(:,end)]; elerr_p=elerror;
error_nodes = [14 16 18 20 13 27 22]; % midpoint face nodes and centroid node
%
% recompute the Q2 matrices
% reconstruct the element matrices
% construct the integration rule
ngpt=3;
[oneg,onew] = gausspoints_oned(ngpt);
[s,t,l,wt] = gausspoints_threed(oneg,onew);
nngpt=ngpt^3; ng=ngpt;
tic
%
% inner loop over elements

xl_v = x(ev); yl_v = y(ev); zl_v = z(ev);
aez = zeros(nel,7,7);
% loop over Gauss points
for igpt = 1:nngpt
    sigpt=s(igpt); tigpt=t(igpt); ligpt=l(igpt); wght=wt(igpt);
    % evaluate derivatives etc
    [~,invjac_v,~,~,~,~] = deriv3D(sigpt,tigpt,ligpt,xl_v,yl_v,zl_v);
    [~,dpsidx_v,dpsidy_v,dpsidz_v] = qderiv3D(sigpt,tigpt,ligpt,xl_v,yl_v,zl_v);
    aez = aez + wght*(dpsidx_v(:,error_nodes).*permute(dpsidx_v(:,error_nodes),[1 3 2]) + dpsidy_v(:,error_nodes).*permute(dpsidy_v(:,error_nodes),[1,3,2])....
        + dpsidz_v(:,error_nodes).*permute(dpsidz_v(:,error_nodes),[1,3,2])).*invjac_v(:);
    % end of Gauss point loop
end

%
% recompute contributions from elements with Dirichlet boundaries
nbdf=length(ebound3D(:,1));
fbdy = zeros(nel,1);
face = zeros(nel,1);
% isolate boundary elements
for el = 1:nbdf
    ee = ebound3D(el,1);
    fbdy(ee) = fbdy(ee)+1; face(ee)=ebound3D(el,2);
end

% three face elements
k3=find(fbdy==3);
nel3b=length(k3);
% loop over two edge elements
for el = 1:nel3b
    el3e=k3(el);
    kk=find(ebound3D(:,1) == el3e);
    faces=ebound3D(kk,2);
    % set up original matrix and RHS vector
    ae=squeeze(aez(el3e,1:7,1:7));
    fe=fez(el3e,:)';
    % set up local coordinates and impose interpolated error as Dirichlet bc
    xl=x(lev(el3e,:)); yl=y(lev(el3e,:)); zl=z(lev(el3e,:));
    [bae,fe] = localbc_p3D(ae,fe,faces,xl,yl,zl);
    % solve local problem
    err=bae\fe;
    elerr_p(el3e,1) = err'*fe;
    
end
% end of element loop
%
%
% two face elements
k2=find(fbdy==2);
nel2b=length(k2);
% loop over two face elements
for el = 1:nel2b
    el2e=k2(el);
    kk=find(ebound3D(:,1) == el2e);
    faces=ebound3D(kk,2);
    % set up original matrix and RHS vector
    ae=squeeze(aez(el2e,1:7,1:7));
    fe=fez(el2e,:)';
    % set up local coordinates and impose interpolated error as Dirichlet bc
    xl=x(lev(el2e,:)); yl=y(lev(el2e,:)); zl=z(lev(el2e,:));
    [bae,fe] = localbc_p3D(ae,fe,faces,xl,yl,zl);
    % solve local problem
    err=bae\fe;
    elerr_p(el2e,1) = err'*fe;
end
% end of element loop
%
% one face elements
k1=find(fbdy==1);
nel1b=length(k1);
% loop over one edge elements
for el = 1:nel1b
    el1e=k1(el);
    kk=find(ebound3D(:,1) == el1e);
    faces=ebound3D(kk,2);
    % set up original matrix and RHS vector
    fe=fez(el1e,:)';
    ae=squeeze(aez(el1e,1:7,1:7));
    % set up local coordinates and impose interpolated error as Dirichlet bc
    xl=x(lev(el1e,:)); yl=y(lev(el1e,:)); zl=z(lev(el1e,:));
    [bae,fe] = localbc_p3D(ae,fe,faces,xl,yl,zl);
    % solve local problem
    err=bae\fe;
    elerr_p(el1e,1) = err'*fe;
end
% end of element loop
%
etime=toc;
fprintf('error boundary correction took %6.3e seconds\n',etime)
return
