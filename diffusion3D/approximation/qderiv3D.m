function [psi,dpsidx,dpsidy,dpsidz] = qderiv3D(s,t,l,xl,yl,zl)
%QDERIV3D evaluates derivatives of triquadratic shape functions
%   [psi,dpsidx,dpsidy,dpsidz] = qderiv3D(s,t,l,xl,yl);
%   input
%          s         reference element x coordinate
%          t         reference element y coordinate
%          l         reference element z coordinate
%          xl        physical element x vertex coordinates
%          yl        physical element y vertex coordinates
%          zl        physical element z vertex coordinates
%   output
%          psi       elementwise shape functions
%          dpsidx    x derivatives of psi
%          dpsidy    y derivatives of psi
%          dpsidz    z derivatives of psi
% IFISS function: GP; 9 June 2022.
% Copyright (c)  2022  G.Papanikos,  C.E. Powell, D.J. Silvester

nel=length(xl(:,1));
zero_v = zeros(nel,1);
one_v = ones(nel,1);
%
% evaluate trilinear shape functions
[phi_e,dphids,dphidt,dphidl] = shape3D(s,t,l);
% evaluate triquadratic shape functions
[psi_e,dpsids,dpsidt,dpsidl] = qshape3D(s,t,l);
% local derivatives
dxds = zero_v;
dxdt = zero_v;
dxdl = zero_v;
dyds = zero_v;
dydt = zero_v;
dydl = zero_v;
dzds = zero_v;
dzdt = zero_v;
dzdl = zero_v;

jac = zero_v;
invjac = zero_v;
%
for ivtx = 1:8
    dxds(:) = dxds(:) + xl(:,ivtx) .* one_v*dphids(ivtx);
    dxdt(:) = dxdt(:) + xl(:,ivtx) .* one_v*dphidt(ivtx);
    dxdl(:) = dxdl(:) + xl(:,ivtx) .* one_v*dphidl(ivtx);
    
    dyds(:) = dyds(:) + yl(:,ivtx) .* one_v*dphids(ivtx);
    dydt(:) = dydt(:) + yl(:,ivtx) .* one_v*dphidt(ivtx);
    dydl(:) = dydl(:) + yl(:,ivtx) .* one_v*dphidl(ivtx);
    
    dzds(:) = dzds(:) + zl(:,ivtx) .* one_v*dphids(ivtx);
    dzdt(:) = dzdt(:) + zl(:,ivtx) .* one_v*dphidt(ivtx);
    dzdl(:) = dzdl(:) + zl(:,ivtx) .* one_v*dphidl(ivtx);
    
end
%
jac(:) =  dxds(:).*(dydt(:).*dzdl(:) - dydl(:).*dzdt(:))...
    -dxdt(:).*(dyds(:).*dzdl(:) - dydl(:).*dzds(:))...
    +dxdl(:).*(dyds(:).*dzdt(:) - dydt(:).*dzds(:));
%
% check element Jacobian
if any(jac < 1e-9)
    fprintf('Bad element warning ...\n')
    if any(jac <= 0.0)
        error('singular Jacobian ... Aborted ...')
    end
end
invjac(:) = one_v ./ jac(:);
%

psi = psi_e.*one_v;

alpha = dydt(:).*dzdl(:) - dydl(:).*dzds(:);
beta =  dxdt(:).*dzdl(:) - dxdl(:).*dzdt(:);
gamma = dxdt(:).*dydl(:) - dxdl(:).*dydt(:);

delta = dyds(:).*dzdl(:) - dydl(:).*dzds(:);
epsilon = dxds(:).*dzdl(:) - dxdl(:).*dzds(:);
zeta = dxds(:).*dydl(:) - dxdl(:).*dyds(:);

eta = dyds(:).*dzdt(:) - dydt(:).*dzds(:);
theta = dxds(:).*dzdt(:) - dxdt(:).*dzds(:);
iota = dxds(:).*dydt(:) - dxdt(:).*dydl(:);

dpsidx = ( dpsids.*one_v.*alpha(:) - dpsidt.*one_v.*beta(:)+dpsidl.*one_v.*gamma(:));
dpsidy= (-dpsids.*one_v.*delta(:) + dpsidt.*one_v.*epsilon(:)-dpsidl.*one_v.*zeta(:));
dpsidz = (dpsids.*one_v.*eta(:) - dpsidt.*one_v.*theta(:)+dpsidl.*one_v.*iota(:));

return