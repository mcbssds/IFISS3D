function EE = energy_norm_error(xyz,ev,FEMsol,sn,qmethod)
% ENERGY_NORM_ERROR computes the energy norm of the error between the
% exact solution and the FEM approximation for the test problems 4,6 with
% exact solution and the H^1 seminorm for the problems 3 and 5
%   EE = energy_norm_error(xyz,ev,q1sol,sn)
%   input
%          xyz        vertex coordinate vector
%          ev         element mapping matrix
%          FEMsol     finite element solution
%          sn         test problems with analytical solution (test problems 3, 4, 5 and 6)
%  outpt
%          EE         energy norm of the error
%
% IFISS function: GP; 9 June 2022.
% Copyright (c)  2022  G.Papanikos,  C.E. Powell, D.J. Silvester

xx = xyz(:,1); yy = xyz(:,2); zz = xyz(:,3);
xl_v = xx(ev);
yl_v = yy(ev);
zl_v = zz(ev);
sl_v = FEMsol(ev);

hx=max(diff(xx)); hy=max(diff(yy));hz=max(diff(zz));

ngpt=3;
[oneg,onew] = gausspoints_oned(ngpt);
[s,t,l,wt] = gausspoints_threed(oneg,onew);
nngpt=ngpt^3; ng=ngpt;


E_numk = zeros(size(ev,1),1);
for igpt = 1:nngpt
    sigpt=s(igpt); tigpt=t(igpt); ligpt=l(igpt);wght=wt(igpt);
    E_numk = E_numk + wght.*eval_error(sigpt,tigpt,ligpt,xl_v,yl_v,zl_v,sl_v,sn,qmethod);
end

EE = sqrt(sum(E_numk));

end

function [jac,xe,ye,ze,soldx ,soldy, soldz] = num_loc_sol(sl_v,s,t,l,xl,yl,zl,qmethod)


[jac,invjac,phi,dphidx,dphidy,dphidz] = deriv3D(s,t,l,xl,yl,zl);
soldx = 0; soldy = 0; soldz = 0;
if qmethod == 1
   
    [nel,~] = size(xl);
    zero_v = zeros(nel,1);
    xe =zero_v; ye =zero_v; ze=zero_v;
    
    for  ivtx=1:8
        soldx = soldx + dphidx(:,ivtx).*sl_v(:,ivtx).*invjac;
        soldy = soldy + dphidy(:,ivtx).*sl_v(:,ivtx).*invjac;
        soldz = soldz + dphidz(:,ivtx).*sl_v(:,ivtx).*invjac;
        
        xe(:) = xe(:) + xl(:,ivtx) .*phi(:,ivtx);
        ye(:) = ye(:) + yl(:,ivtx) .*phi(:,ivtx);
        ze(:) = ze(:) + zl(:,ivtx) .*phi(:,ivtx);
    end
    
else
    [psi,dpsidx,dpsidy,dpsidz] = qderiv3D(s,t,l,xl,yl,zl);
    [nel,~] = size(xl);
    zero_v = zeros(nel,1);
    xe =zero_v; ye =zero_v; ze=zero_v;
    for  ivtx=1:27
        soldx = soldx + dpsidx(:,ivtx).*sl_v(:,ivtx).*invjac;
        soldy = soldy + dpsidy(:,ivtx).*sl_v(:,ivtx).*invjac;
        soldz = soldz + dpsidz(:,ivtx).*sl_v(:,ivtx).*invjac;
        xe(:) = xe(:) + xl(:,ivtx) .*psi(:,ivtx);
        ye(:) = ye(:) + yl(:,ivtx) .*psi(:,ivtx);
        ze(:) = ze(:) + zl(:,ivtx) .*psi(:,ivtx);
    end
    
    
end

end

function Ek = eval_error(s,t,l,x,y,z,sl_v,sn,qmethod)

[jac,xe,ye,ze,soldx ,soldy, soldz] = num_loc_sol(sl_v,s,t,l,x,y,z,qmethod);

if sn  == 11 % derivatives of trilinear exact solution u = -(x-1)(y-1)(z-1)
    ux = -(ye - 1).*(ze - 1);
    uy = -(xe - 1).*(ze - 1);
    uz = -(xe - 1).*(ye - 1);
elseif sn == 12 % derivatives of triquadratic exact solution u = -(x.^2-1)(y.^2-1)(z.^2-1)
    ux = -2.*xe.*(ye.^2 - 1).*(ze.^2 - 1);
    uy = -2.*ye.*(xe.^2 - 1).*(ze.^2 - 1);
    uz = -2.*ze.*(xe.^2 - 1).*(ye.^2 - 1);
elseif sn == 13 % derivatives of tricubic exact solution u = -(x.^3-1)(y.^3-1)(z.^3-1)
    ux = -3.*xe.^2.*(ye.^3 - 1).*(ze.^3 - 1);
    uy = -3.*ye.^2.*(xe.^3 - 1).*(ze.^3 - 1);
    uz = -3.*ze.^2.*(xe.^3 - 1).*(ye.^3 - 1);
elseif sn == 14 % derivatives of triquartic exact solution u = -(x.^4-1)(y.^4-1)(z.^4-1)
    ux = -4.*xe.^3.*(ye.^4 - 1).*(ze.^4 - 1);
    uy = -4.*ye.^3.*(xe.^4 - 1).*(ze.^4 - 1);
    uz = -4.*ze.^3.*(xe.^4 - 1).*(ye.^4 - 1);
else
    error('Wrong test problem!')
end



ex = (ux - soldx).^2;
ey = (uy - soldy).^2;
ez = (uz - soldz).^2;


Ek = (ex+ey+ez).*jac;

end
