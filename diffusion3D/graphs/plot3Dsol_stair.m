function plot3Dsol_stair(sol3D,xyzleft,xright,x,y,z,fig)
%PLOT3DSOL_STAIR plots solution on stair domain
%   plot3Dsol_stair(sol3D,xyzleft,xright,x,y,z,fig);
%   inputs:
%          sol3D           nodal solution vector
%          xyzleft,xright  vertex coordinates
%          x               vector of x-axis interpolation points
%          y               vector of y-axis interpolation points
%          z               vector of z-axis interpolation points
%          fig             figure number
%
% IFISS function: GP; 9 June 2022.
% Copyright (c)  2022  G.Papanikos,  C.E. Powell, D.J. Silvester

fprintf('plotting solution... \n')

[X,Y,Z]=meshgrid(x,y,z);
xyz = [xyzleft;xright];

sol3D = griddata(xyz(:,1),xyz(:,2),xyz(:,3),sol3D,X,Y,Z);

for i=1:size(sol3D,3)
    [II,JJ] = find(X(:,:,i)<0 & Z(:,:,i)<0);
    sol3D(II,JJ,i) =nan;
end

sol3D = permute(sol3D,[3 2 1]);

% create axes and set up initial properties
figure(fig)
subplot(121),contour(X(:,:,round(size(X,3)/2)),Y(:,:,round(size(Y,3)/2)),sol3D(:,:,round(size(Z,3)/2)),20),axis('square')
axis('off'), squarex, title('Finite Element Solution','FontSize',12)

subplot(122),
daspect([1 1 1]);
axis([-1 1 -1 1 -1 1])
xlabel('x') % x-axis label
ylabel('y') % y-axis label
zlabel('z') % z-axis label

hold on

hSlice = slice(x,y,z,sol3D,1/2,1/2,0);
set(hSlice,'EdgeColor','none','FaceColor','interp');
% adjust lighting
camlight;
camlight(-90,0);
lighting gouraud
colormap jet; colorbar;
view(330,30)

return
