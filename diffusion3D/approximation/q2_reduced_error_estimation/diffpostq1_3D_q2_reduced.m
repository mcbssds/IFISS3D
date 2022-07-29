function [errorsq_ele,xx,fe,ae] = diffpostq1_3D_q2_reduced(xyz,ev,ebound3D,q1sol3D,fcx,hx,hy,hz)
%DIFFPOSTQ1_3D_Q2_reduced local Poisson error estimator for Q1 solution
%using only Q2 base functions in the midpoint of the face and centroid point
%   [errorsq_ele,elerr,fe,ae] = diffpostq1_q2_reduced(xyz,ev,ebound3D,q1sol3D,fcx,hx,hy);
%   input
%          xyz           vertex coordinate vector
%          ev            element mapping matrix
%          ebound3D      element edge boundary matrix
%          q1sol3D       Q1 solution vector
%          fcx           element face connectivity array
%          hx,hy,hz      element mesh sizes
%   output
%          errorsq_ele  element error estimate
%          elerr        elementwise error estimate
%          fe           elementwise rhs vectors
%          ae           LDLT factorized element matrices
%
%   calls functions q1fluxjmps3D, gausspoints_oned, gausspoints_threed,
%   deriv3D, qderiv3D gauss_source3D
%   IFISS function: GP; 09 June 2022
% Copyright (c)  2022 G.Papanikos, C.E. Powell, D.J. Silvester

fprintf('computing Q1 error estimator...  \n')
x=xyz(:,1); y=xyz(:,2); z=xyz(:,3);
nel=length(ev(:,1));
errorsq_ele = zeros(nel,1);
%
% construct the integration rule
ngpt=3;
[oneg,onew] = gausspoints_oned(ngpt);
[s,t,l,wt] = gausspoints_threed(oneg,onew);
nngpt=ngpt^3; ng=ngpt;
tic

error_nodes = [14 16 18 20 13 27 22]; % midpoint face nodes and centroid node
%
% inner loop over elements

xl_v = x(ev);  yl_v = y(ev); zl_v = z(ev);
ae = zeros(nel,7,7); elerr=zeros(7,nel);
fe = zeros(nel,7);
% loop over Gauss points
for igpt = 1:nngpt
    sigpt=s(igpt);
    tigpt=t(igpt);
    ligpt= l(igpt);
    wght=wt(igpt);
    % evaluate derivatives etc
    [jac_v,invjac_v,~,~,~,~] = deriv3D(sigpt,tigpt,ligpt,xl_v,yl_v,zl_v);
    rhs_v = gauss_source3D(sigpt,tigpt,ligpt,xl_v,yl_v,zl_v);
    [psi_v,dpsidx_v,dpsidy_v,dpsidz_v] = qderiv3D(sigpt,tigpt,ligpt,xl_v,yl_v,zl_v);

    fe = fe + wght*rhs_v(:).*psi_v(:,error_nodes).*jac_v(:);
    ae = ae + wght*(dpsidx_v(:,error_nodes).*permute(dpsidx_v(:,error_nodes),[1 3 2]) + dpsidy_v(:,error_nodes).*permute(dpsidy_v(:,error_nodes),[1,3,2])....
        + dpsidz_v(:,error_nodes).*permute(dpsidz_v(:,error_nodes),[1,3,2])).*invjac_v(:);
    
    % end of Gauss point loop
end
%
% include face jumps (evaluated at mid nodes on each face)
njmp_midpoint  = q1fluxjmps3D(q1sol3D,fcx,xyz,ev,ebound3D,0);
%----------------------------------

fe(:,1) = fe(:,1)- 2.*njmp_midpoint(:,1).*hx(:).*hy(:)./9;  % midnode on face 1
fe(:,2) = fe(:,2)- 2.*njmp_midpoint(:,2).*hy(:).*hz(:)./9;  % midnode on face 2
fe(:,3) = fe(:,3)- 2.*njmp_midpoint(:,3).*hy(:).*hx(:)./9;  % midnode on face 3
fe(:,4) = fe(:,4)- 2.*njmp_midpoint(:,4).*hy(:).*hz(:)./9;  % midnode on face 4
fe(:,5) = fe(:,5)- 2.*njmp_midpoint(:,5).*hx(:).*hz(:)./9;  % midnode on face 5
fe(:,6) = fe(:,6) -2.*njmp_midpoint(:,6).*hx(:).*hz(:)./9;  % midnode on face 6


% LDLT factorization
nn=7;
dd=zeros(nel,nn); rr=zeros(nel,nn);
for kk=1:nn-1
    for pp=1:kk-1;
        rr(1:nel,pp)=dd(1:nel,pp).*ae(1:nel,kk,pp);
    end
    dd(1:nel,kk)=ae(1:nel,kk,kk);
    for pp=1:kk-1;
        dd(1:nel,kk)= dd(1:nel,kk)- ae(1:nel,kk,pp).*rr(1:nel,pp);
    end
    for ii=kk+1:nn
        for pp=1:kk-1;
            ae(1:nel,ii,kk)=ae(1:nel,ii,kk)-ae(1:nel,ii,pp).*rr(1:nel,pp);
        end
        ae(1:nel,ii,kk)=ae(1:nel,ii,kk)./dd(1:nel,kk);
    end
end
for pp=1:nn-1;
    rr(1:nel,pp)=dd(1:nel,pp).*ae(1:nel,nn,pp);
end
dd(1:nel,nn)=ae(1:nel,nn,nn);
for pp=1:nn-1;
    dd(1:nel,nn)= dd(1:nel,nn)- ae(1:nel,nn,pp).*rr(1:nel,pp);
end
% overwrite diagonal entries
for kk=1:nn
    ae(1:nel,kk,kk)= dd(1:nel,kk);
end
%  forward-backward substitutions ...
xx = element_lusolve(ae,fe);
elerr=xx';

%
for ivtx=1:7,
    errorsq_ele(:) = errorsq_ele(:) + fe(:,ivtx) .* elerr(ivtx,:)';
end
etime=toc;
fprintf('error estimation took %6.3e seconds\n',etime)
return