function [A,Q,f] = femq1_diff3D(xyz,ev)
%FEMQ1_DIFF3D vectorized trilinear coefficient matrix generator
%   [A,Q,f] = femq1_diff3D(xy,ev);
%   inputs:
%          xyz        vertex coordinates
%          ev         element mapping matrix
%   outputs:
%          A          stiffness matrix
%          Q          mass matrix
%          f          rhs vector
%
% Natural boundary conditions applied. Dirichlet conditions
% must be explicitly enforced by calling function nonzerobc3D.
% IFISS function: GP; 9 June 2022.
% Copyright (c) 2021 G. Papanikos, C.E. Powell, D.J. Silvester

x=xyz(:,1); y=xyz(:,2); z=xyz(:,3);
nvtx=length(x);
[nel,nvt]=size(ev);
lx=max(x)-min(x); ly=max(y)-min(y);lz=max(z)-min(z);
hx=max(diff(x)); hy=max(diff(y));hz=max(diff(z));
fprintf('setting up Q1 diffusion matrices...  ')
%
% initialise global matrices
A = sparse(nvtx,nvtx);
Q = sparse(nvtx,nvtx);
f = zeros(nvtx,1);
%

% construct the integration rule
ngpt=3;                             % 2x2x2 Gauss points %(s,t,l)_i = ((-1)^k gpt,(-1)^j gpt,(-1)^r gpt), i = 4(k-1) +2(j-1) +r , k,j,r =1,2
[oneg,onew] = gausspoints_oned(ngpt);
[s,t,l,wt] = gausspoints_threed(oneg,onew);
nngpt=ngpt^3; ng=ngpt;

% inner loop over elements

xl_v = x(ev); yl_v = y(ev); zl_v = z(ev);

ae = zeros(nel,8,8);
re = zeros(nel,8,8);
fe = zeros(nel,8);

%  loop over Gauss points
for igpt = 1:nngpt
    sigpt=s(igpt);  tigpt=t(igpt); ligpt=l(igpt);  wght =wt(igpt);
    %  evaluate derivatives etc
    [jac,invjac,phi,dphidx,dphidy,dphidz] = deriv3D(sigpt,tigpt,ligpt,xl_v,yl_v,zl_v);
    
    rhs = gauss_source3D(sigpt,tigpt,ligpt,xl_v,yl_v,zl_v);
    fe = fe + wght*rhs(:).*phi.*jac(:);
    ae = ae + wght*(dphidx.*permute(dphidx,[1 3 2]) + dphidy.*permute(dphidy,[1,3,2]) + dphidz.*permute(dphidz,[1,3,2]) ).*invjac(:);
    re = re + wght*phi.*permute(phi,[1 3 2]);
    % end of Gauss point loop
end

% perform assembly of global matrix and source vector 
for krow=1:8
    nrow=ev(:,krow);
    for kcol=1:8
        ncol=ev(:,kcol);
        A = A + sparse(nrow,ncol,ae(:,krow,kcol),nvtx,nvtx);
        Q = Q + sparse(nrow,ncol,re(:,krow,kcol),nvtx,nvtx);
    end
    f(nrow,1) = f(nrow,1) + fe(:,krow);
end

fprintf('done\n')
return
