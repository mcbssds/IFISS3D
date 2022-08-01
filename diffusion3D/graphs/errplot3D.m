function errplot3D(sol3D,eldata3D,ev,xyz,x,y,z,fig)
%ERRPLOT3D plots solution and error on cube domain
%   errplot3D(sol3D,eldata3D,ev,xyz,x,y,z,fig);
%   inputs:
%          sol3D        nodal solution vector
%          eldata3D     element error vector
%          ev           element mapping matrix
%          xyz          vertex coordinates
%          x            vector of x-axis interpolation points
%          y            vector of y-axis interpolation points
%          z            vector of z-axis interpolation points
%          fig          figure number
%
% IFISS function: GP; 9 June 2022.
% Copyright (c)  2022  G.Papanikos,  C.E. Powell, D.J. Silvester

fprintf('plotting solution and estimated errors... ')

[X,Y,Z]=meshgrid(x,y,z);

sol3D = reshape(sol3D,[size(x,1),size(y,1),size(z,1)]);
xmin = min(x(:)); xmax =  max(x(:)); xmean = mean(x(:)); 
ymin = min(y(:)); ymax =  max(y(:)); ymean =  mean(y(:));
zmin = min(z(:)); zmax =  max(z(:)); zmean =  mean(z(:));
% create axes and set up initial properties
figure(fig)
subplot(221),contour(X(:,:,round(size(X,3)/2)),Y(:,:,round(size(Y,3)/2)),sol3D(:,:,round(size(Z,3)/2)),20),axis('square')
axis('off'), squarex, title('Finite Element Solution','FontSize',12)
axis('square')
subplot(222),
daspect([1 1 1]);
axis([xmin xmax ymin ymax zmin zmax])
xlabel('x') % x-axis label
ylabel('y') % y-axis label
zlabel('z') % z-axis label
hold on
hSlice = slice(x,y,z,sol3D,xmean,ymean,zmin);
set(hSlice,'EdgeColor','none','FaceColor','interp');

hSlice = slice(x,y,z,sol3D, (xmax - xmean)/2,(ymax - ymean)/2,zmean);
set(hSlice,'EdgeColor','none','FaceColor','interp');

% adjust lighting
camlight;
camlight(-90,0);
lighting gouraud
colormap jet; colorbar;
view(330,30)

xx=xyz(:,1); yy=xyz(:,2); zz=xyz(:,3);
nel=length(eldata3D);

xl = xx(ev);
yl = yy(ev);
zl = zz(ev);
xc(:,1) = 0.125*sum(xl,2);
xc(:,2) = 0.125*sum(yl,2);
zc(:,3) = 0.125*sum(zl,2);

% interpolate to a cartesian product mesh
x=0.5*(x(1:end-1)+x(2:end));
y=0.5*(y(1:end-1)+y(2:end));
z=0.5*(z(1:end-1)+z(2:end));
[X,Y,Z]=meshgrid(x,y,z);
xyzsol = griddata(xc(:,1),xc(:,2),zc(:,3),eldata3D,X,Y,Z);

subplot(223),contour(X(:,:,round(size(X,3)/2)),Y(:,:,round(size(Y,3)/2)),xyzsol(:,:,round(size(Z,3)/2)),15),axis('square')
title('Estimated Error','FontSize',12)
subplot(224)
daspect([1 1 1]);
axis([xmin xmax ymin ymax zmin zmax])
xlabel('x') % x-axis label
ylabel('y') % y-axis label
zlabel('z') % z-axis label
hold on
hSlice = slice(x,y,z,xyzsol,xmean,ymean,zmin);
set(hSlice,'EdgeColor','none','FaceColor','interp');

hSlice = slice(x,y,z,xyzsol,(xmax - xmean)/2,(ymax - ymean)/2,zmean);
set(hSlice,'EdgeColor','none','FaceColor','interp');

% adjust lighting
camlight;
camlight(-90,0);
lighting gouraud
colormap jet; colorbar;
view(330,30)
fprintf('done\n')

return
