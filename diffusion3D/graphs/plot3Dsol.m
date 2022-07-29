function plot3Dsol(sol3D,x,y,z,fig)
%PLOT3Dsol plots the 3D solution and error on cube domain
%   plot3DQ2(sol3D,x,y,z,fig);
%   input
%          sol3D        nodal solution vector 
 
%          x            vector of x-axis interpolation points
%          y            vector of y-axis interpolation points
%          z            vector of z-axis interpolation points
%          fig          figure number
%
% IFISS function: GP; 9 June 2022.
% Copyright (c)  2022  G.Papanikos,  C.E. Powell, D.J. Silvester
fprintf('plotting solution...\n')

[X,Y,Z]=meshgrid(x,y,z);
xmin = min(x(:)); xmax =  max(x(:)); xmean = mean(x(:)); 
ymin = min(y(:)); ymax =  max(y(:)); ymean =  mean(y(:));
zmin = min(z(:)); zmax =  max(z(:)); zmean =  mean(z(:));
sol3D = reshape(sol3D,[size(x,1),size(y,1),size(z,1)]);


%% create an axes and setup initial properties
figure(fig)
subplot(121),contour(X(:,:,round(size(X,3)/2)),Y(:,:,round(size(Y,3)/2)),sol3D(:,:,round(size(Z,3)/2)),20),axis('square')
axis('off'), squarex, title('Finite Element Solution','FontSize',12)

subplot(122),
daspect([1 1 1]);
axis([xmin xmax ymin ymax zmin zmax])
xlabel('x') % x-axis label
ylabel('y') % y-axis label
zlabel('z') % z-axis label

hold on
hSlice = slice(x,y,z,sol3D,xmean,ymean,zmin);
set(hSlice,'EdgeColor','none','FaceColor','interp');

hSlice = slice(x,y,z,sol3D,(xmax - xmean)/2,(ymax - ymean)/2,zmean);
set(hSlice,'EdgeColor','none','FaceColor','interp');
%% adjust lighting
camlight;
camlight(-90,0);
lighting gouraud
colorbar
view(330,30)

return