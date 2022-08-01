function [x,y,z,xyz] = q2grid3D(x,y,z,xyz,mv,bound);
%Q2GRID3D triquadratic element grid generator
%   [x,y,z,xyz] = q2grid3D(x,y,z,xyz,mv,bound);
%   inputs:
%          x          x coordinate vector
%          y          y coordinate vector
%          z          z coordinate vector
%          xyz        vertex coordinates
%          mv         Q2 macroelement mapping matrix
%          bound      boundary vertex vector
%   outputs:
%
% postpocesses Q2 element partitioning information to
% give standard approximation in the case of stretched grids
% IFISS function: GP; 9 June 2022; .
% Copyright (c)  G.Papanikos, C.E. Powell, D.J. Silvester

xx=xyz(:,1); yy=xyz(:,2); zz = xyz(:,3); nvtx=length(xx);
%
% recompute mid-side points in the case of stretched grids
% y-direction
yv=yy; ny=length(y);
for k=2:2:ny;
    yold=y(k); ynew=0.5*(y(k+1)+y(k-1));
    l=find(yy==yold); yv(l)=ynew; y(k)=ynew;
end
% x-direction
xv=xx; nx=length(x);
for k=2:2:nx;
    xold=x(k); xnew=0.5*(x(k+1)+x(k-1));
    l=find(xx==xold); xv(l)=xnew; x(k)=xnew;
end

% z-direction
zv=zz; nx=length(z);
for k=2:2:nx;
    zold=z(k); znew=0.5*(z(k+1)+z(k-1));
    l=find(zz==zold); zv(l)=znew; z(k)=znew;
end


xyz=[xv,yv,zv];

return
