function [errorsq_ele,xx,fe,ae] = diffpostq1_q2_3D(xyz,ev,ebound3D,q1sol3D,fcx,hx,hy,hz)
%DIFFPOSTQ1_Q2_3D local Poisson error estimator for Q1 solution using full Q2 bubble functions
%  [errorsq_ele,xx,fe,ae] = diffpostq1_q2_3D(xyz,ev,ebound3D,q1sol3D,fcx,hx,hy,hz);
%  inputs:
%          xyz           vertex coordinate vector
%          ev            element mapping matrix
%          ebound3D      element edge boundary matrix
%          q1sol3D       Q1 solution vector
%          fcx           element face connectivity array
%          hx,hy,hz      element mesh sizes
%   outputs:
%          errorsq_ele  element error estimate
%          xx           elementwise error estimate
%          fe           elementwise rhs vectors
%          ae           LDLT factorized element matrices
%
%   calls functions q1fluxjmps3D, q1fluxjmps3Dv2, gausspoints_oned, gausspoints_threed
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

xl_v = x(ev); yl_v = y(ev); zl_v = z(ev);

ae = zeros(nel,19,19); elerr=zeros(19,nel);
fe = zeros(nel,19);
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

    fe = fe + wght*rhs_v(:).*psi_v(:,9:27).*jac_v(:);
    ae = ae + wght*(dpsidx_v(:,9:27).*permute(dpsidx_v(:,9:27),[1 3 2]) + dpsidy_v(:,9:27).*permute(dpsidy_v(:,9:27),[1,3,2])....
        + dpsidz_v(:,9:27).*permute(dpsidz_v(:,9:27),[1,3,2])).*invjac_v(:);
    % end of Gauss point loop
end
%
% include face jumps (evaluated at the nodes on each face)
njmp_midpoint  = q1fluxjmps3D(q1sol3D,fcx,xyz,ev,ebound3D,0);  % find flux jump on the mid face node
[njmp1,njmp2,njmp3,njmp4]  = q1fluxjmps3Dv2(q1sol3D,fcx,xyz,ev,ebound3D,0,0,0);  % find flux jumps on face edge nodes
%----------------------------------------------------------
%------------------node numbering -------------------------
% 9 10 11 12 13  14  15 16 17 18 19 20 21 22 23 24 25 26 27
% 1  2  3  4  5   6   7  8  9 10 11 12 13 14 15 16 17 18 19

% face 1
fe(:,1)  = fe(:,1)   - njmp1(:,1).*hx(:).*hy(:)./18;
fe(:,7)  = fe(:,7)   - njmp2(:,1).*hx(:).*hy(:)./18;
fe(:,15) = fe(:,15)  - njmp3(:,1).*hx(:).*hy(:)./18;
fe(:,13) = fe(:,13)  - njmp4(:,1).*hx(:).*hy(:)./18;
fe(:,6)  = fe(:,6)   - 2.*njmp_midpoint(:,1).*hx(:).*hy(:)./9; %midnode

%face 2
fe(:,7)  = fe(:,7)   - njmp1(:,2).*hy(:).*hz(:)./18;
fe(:,16) = fe(:,16)  - njmp2(:,2).*hy(:).*hz(:)./18;
fe(:,9)  = fe(:,9)   - njmp3(:,2).*hy(:).*hz(:)./18;
fe(:,2)  = fe(:,2)   - njmp4(:,2).*hy(:).*hz(:)./18;
fe(:,8)  = fe(:,8)   - 2.*njmp_midpoint(:,2).*hy(:).*hz(:)./9; % midnode

%face 3
fe(:,3)  = fe(:,3)   - njmp1(:,3).*hx(:).*hy(:)./18;
fe(:,9)  = fe(:,9)   - njmp2(:,3).*hx(:).*hy(:)./18;
fe(:,17) = fe(:,17)  - njmp3(:,3).*hx(:).*hy(:)./18;
fe(:,11) = fe(:,11)  - njmp4(:,3).*hx(:).*hy(:)./18;
fe(:,10) = fe(:,10)  - 2.*njmp_midpoint(:,3).*hy(:).*hx(:)./9; % midnode

%face 4
fe(:,13) = fe(:,13)  - njmp1(:,4).*hy(:).*hz(:)./18;
fe(:,18) = fe(:,18) -  njmp2(:,4).*hy(:).*hz(:)./18;
fe(:,11) = fe(:,11)  - njmp3(:,4).*hy(:).*hz(:)./18;
fe(:,4)  = fe(:,4)   - njmp4(:,4).*hy(:).*hz(:)./18;
fe(:,12) = fe(:,12)  - 2.*njmp_midpoint(:,4).*hy(:).*hz(:)./9; % midnode

%face 5
fe(:,1)  = fe(:,1)  - njmp1(:,5).*hx(:).*hz(:)./18;
fe(:,2)  = fe(:,2)  - njmp2(:,5).*hx(:).*hz(:)./18;
fe(:,3)  = fe(:,3)  - njmp3(:,5).*hx(:).*hz(:)./18;
fe(:,4)  = fe(:,4)  - njmp4(:,5).*hx(:).*hz(:)./18;
fe(:,5)  = fe(:,5)  - 2.*njmp_midpoint(:,5).*hx(:).*hz(:)./9; % midnode


%face 6
fe(:,15) = fe(:,15) - njmp1(:,6).*hx(:).*hz(:)./18;
fe(:,16) = fe(:,16) - njmp2(:,6).*hx(:).*hz(:)./18;
fe(:,17) = fe(:,17) - njmp3(:,6).*hx(:).*hz(:)./18;
fe(:,18) = fe(:,18) - njmp4(:,6).*hx(:).*hz(:)./18;
fe(:,19) = fe(:,19) - 2.*njmp_midpoint(:,6).*hx(:).*hz(:)./9; % midnode


% vectorized code
% LDLT factorization
nn=19;
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
for ivtx=1:19,
    errorsq_ele(:) = errorsq_ele(:) + fe(:,ivtx) .* elerr(ivtx,:)';
end
etime=toc;
fprintf('error estimation took %6.3e seconds\n',etime)
return
