function [A,Q,f] = femq2_diff3D(xyz,mv)
%FEMQ2_DIFF vectorized triquadratic coefficient matrix generator
%   [A,Q,f] = femq2_diff3D(xyz,mv);
%   input
%          xyz         nodal coordinate vector
%          mv         element mapping matrix
%   output
%          A          stiffness matrix
%          Q          mass matrix
%          f          rhs vector
%
%   Natural boundary conditions apply. Dirichlet conditions
%   must be explicitly enforced by calling function nonzerobc.
%   IFISS function: GP; 9 June 2022.
% Copyright (c) 2021 G. Papanikos, C. E Powell, D.J. Silvester

x=xyz(:,1); y=xyz(:,2); z= xyz(:,3);
nvtx=length(x);
nel=length(mv(:,1));
lx=max(x)-min(x); ly=max(y)-min(y);lz=max(z)-min(z);
hx=max(diff(x)); hy=max(diff(y)); hz =max(diff(z));
fprintf('setting up Q2 diffusion matrices...  ')
%
% initialise global matrices
A = sparse(nvtx,nvtx);
Q = sparse(nvtx,nvtx);
f = zeros(nvtx,1);
%
% construct the integration rule
ngpt=3;
[oneg,onew] = gausspoints_oned(ngpt);
[s,t,l,wt] = gausspoints_threed(oneg,onew);
nngpt=ngpt^3; ng=ngpt;
%
xl_v = x(mv(:,1:8));  yl_v = y(mv(:,1:8)); zl_v = z(mv(:,1:8));

ae = zeros(nel,27,27);
re = zeros(nel,27,27);
fe = zeros(nel,27);
% loop over Gauss points
for igpt = 1:nngpt
    sigpt=s(igpt);
    tigpt=t(igpt);
    ligpt=l(igpt);
    wght=wt(igpt);
    %  evaluate derivatives etc
    [jac,invjac,phi,dphidx,dphidy,dphidz] = deriv3D(sigpt,tigpt,ligpt,xl_v,yl_v,zl_v);
    rhs = gauss_source3D(sigpt,tigpt,ligpt,xl_v,yl_v,zl_v);
    [psi,dpsidx,dpsidy,dpsidz] = qderiv3D(sigpt,tigpt,ligpt,xl_v,yl_v,zl_v);
    
    fe = fe + wght*rhs(:).*psi.*jac(:);
    ae = ae + wght*(dpsidx.*permute(dpsidx,[1 3 2]) + dpsidy.*permute(dpsidy,[1,3,2]) + dpsidz.*permute(dpsidz,[1,3,2]) ).*invjac(:);
    re = re + wght*psi.*permute(psi,[1 3 2]);
    % end of Gauss point loop
end
%
% perform assembly of global matrix  and source vector
for krow=1:27
    nrow=mv(:,krow);
    for kcol=1:27
        ncol=mv(:,kcol);
        A = A + sparse(nrow,ncol,ae(:,krow,kcol),nvtx,nvtx);
        Q = Q + sparse(nrow,ncol,re(:,krow,kcol),nvtx,nvtx);
    end
    f(nrow,1) = f(nrow,1) + fe(:,krow);
end

fprintf('done\n')
return

