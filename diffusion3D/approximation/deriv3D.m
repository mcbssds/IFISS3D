function [jac,invjac,phi,dphidx,dphidy,dphidz] = deriv3D(s,t,l,xl,yl,zl)
%DERIV3D evaluates derivatives of trilinear shape functions
%   [jac,invjac,phi,dphidx,dphidy,dphidz] = deriv(s,t,l,xl,yl,zl);
%   input
%          s         reference element x coordinate
%          t         reference element y coordinate
%          l         reference element z coordinate

%          xl        physical element x vertex coordinates
%          yl        physical element y vertex coordinates
%          zl        physical element z vertex coordinates
%   output
%          jac       elementwise jacobian (evaluated at (s,t,l))
%          invjac    elementwise inverse of jacobian
%          phi       elementwise shape functions
%          dphidx    x derivatives of phi
%          dphidy    y derivatives of phi
%          dphidz    z derivatives of phi
% IFISS function: GP; 9 June 2022.
% Copyright (c)  2022  G.Papanikos,  C.E. Powell, D.J. Silvester

ngpt = size(s,1);

[nel,nv] = size(xl);
zero_v = zeros(nel,ngpt);
one_v = ones(nel,ngpt);
xxl = repmat(xl,1,ngpt);
yyl = repmat(yl,1,ngpt);
zzl = repmat(zl,1,ngpt);
one_vv = ones(nel,nv*ngpt);

% evaluate shape functions
[phi_e,dphids,dphidt,dphidl] = shape3D(s,t,l);
%

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

for ivtx = 1:8
    phi(:,ivtx) = phi_e(ivtx)*one_v;
    
    alpha = dydt(:).*dzdl(:) - dydl(:).*dzds(:);
    beta = dxdt(:).*dzdl(:) - dxdl(:).*dzdt(:);
    gamma = dxdt(:).*dydl(:) - dxdl(:).*dydt(:);
    
    delta = dyds(:).*dzdl(:) - dydl(:).*dzds(:);
    epsilon = dxds(:).*dzdl(:) - dxdl(:).*dzds(:);
    zeta = dxds(:).*dydl(:) - dxdl(:).*dyds(:);
    
    eta = dyds(:).*dzdt(:) - dydt(:).*dzds(:);
    theta = dxds(:).*dzdt(:) - dxdt(:).*dzds(:);
    iota = dxds(:).*dydt(:) - dxdt(:).*dydl(:);
    
    dphidx(:,ivtx) = ( dphids(ivtx).*alpha(:) - dphidt(ivtx).*beta(:)+dphidl(ivtx).*gamma(:));
    dphidy(:,ivtx) = (-dphids(ivtx).*delta(:) + dphidt(ivtx).*epsilon(:)-dphidl(ivtx).*zeta(:));
    dphidz(:,ivtx) = (dphids(ivtx).*eta(:) - dphidt(ivtx).*theta(:)+dphidl(ivtx).*iota(:));
end
return