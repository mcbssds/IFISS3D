function plot3Dsol_borehole(sol3D,xyz,x,y,z,bd,fig)
%PLOTSOL3D_BOREHOLE plots solution on borehole domain
%   plotsol3D_borehole(sol3D,x,y,z,bd,fig);
%   inputs:
%          sol3D        nodal solution vector 
%          eldata3D     element error vector
%          ev           element mapping matrix
%          xyz          vertex coordinate vector  
%          x            vector of x-axis interpolation points
%          y            vector of y-axis interpolation points
%          z            vector of z-axis interpolation points
%          bd           borehole semi-width
%          fig          figure number
%
% IFISS function: GP; 9 June 2022.
% Copyright (c)  2022  G.Papanikos,  C.E. Powell, D.J. Silvester

fprintf('plotting solution... ')

[X,Y,Z]=meshgrid(x,y,z);
sol3D = griddata(xyz(:,1),xyz(:,2),xyz(:,3),sol3D,Y,Z,X);

xmin = min(x(:)); xmax =  max(x(:)); xmean = mean(x(:)); 
ymin = min(y(:)); ymax =  max(y(:)); ymean =  mean(y(:));
zmin = min(z(:)); zmax =  max(z(:)); zmean =  mean(z(:));

for i=1:size(sol3D,3)
    [II,JJ] = find(X(:,:,i)>=-bd & X(:,:,i)<=bd & Y(:,:,i)>=-bd & Y(:,:,i)<=bd & Z(:,:,i)>0 & Z(:,:,i)<=1);
    sol3D(II,JJ,i) =nan; 
end


%% create axes and setup initial properties
figure(fig)
subplot(121),contour(X(:,:,round(size(X,3)/2)),Y(:,:,round(size(Y,3)/2)),sol3D(:,:,round(size(Z,3)/2)),20),axis('square')
axis('off'), squarex, title('Finite Element Solution','FontSize',12)
axis('square')
subplot(122),
daspect([1 1 1]);
axis([xmin xmax ymin ymax zmin zmax])
xlabel('x') % x-axis label
ylabel('z') % z-axis label
zlabel('y') % y-axis label
hold on
hSlice = slice(x,y,z,sol3D,xmean,ymean,zmean);
set(hSlice,'EdgeColor','none','FaceColor','interp');

hSlice = slice(x,y,z,sol3D, (xmax - xmean)/2,(ymax - ymean)/2,zmean);
set(hSlice,'EdgeColor','none','FaceColor','interp');
% adjust lighting
camlight;
camlight(-90,0);
lighting gouraud
colormap jet; colorbar;
view(330,30)
fprintf('done\n')
return
