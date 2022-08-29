function cube_domain(Dmax,nc)
%CUBE_DOMAIN cube domain Q2-element grid generator
%    cube_domain(Dmax,nc);
%    inputs:
%    Dmax    limits of domain   [-Dmax, Dmax] 
%    nc      grid parameter    2^nc x 2^nc x 2^nc cube
%
% grid defining data is saved to the file: cube_grid.mat
% IFISS function: DJS; 29 July 2022
% Copyright (c) G. Papanikos, C.E. Powell, D.J. Silvester

if nc<2,
    error('illegal parameter choice, try again.'),
end
if Dmax<0.1 || Dmax ==0
error('domain is too small, try again.')
end

grid_type=default('Uniform/stretched grid (1/2) (default is uniform)',1);
n=2^nc;
np=n/(Dmax*2);

% compute (x,y,z) coordinates of vertices
% y-direction
if grid_type==2
    hmax=Dmax*nc/(2^(nc+1));
    
    x1=-Dmax;
    x2=-2*hmax;
    x3=2*hmax;
    x4=Dmax;
    
    nx1=2^(nc-1)-1;
    nx2=2;
    nx3=2^(nc-1)-1;
    
    y1=-Dmax;
    y2=-2*hmax;
    y3=2*hmax;
    y4=Dmax;
    ny1=2^(nc-1)-1;
    ny2=2;
    ny3=2^(nc-1)-1;
    
    y=subint(y1,y2,y3,y4,ny1,ny2,ny3);
    
    stretch=(y(3)-y(2))/(y(2)-y(1));
    left=-Dmax;
    x=y;
    z=y;
else
    yy=[1/np:1/np:Dmax];
    ypos=[0,yy];
    yneg=-yy(length(yy):-1:1);
    y=[yneg,ypos]';
    left=-Dmax;
    x=y;
    z=y;
    
end

nvtx=(n+1)*(n+1)*(n+1);
[X,Y,Z]=meshgrid(x,y,z);

xx=reshape(Y,nvtx,1);
yy=reshape(Z,nvtx,1);
zz=reshape(X,nvtx,1);

xyz=[xx(:),yy(:),zz(:)];
nn = length(1:2:np*Dmax+np);
v  = 1:2:np*Dmax+np;
xn = 1:np*np*np*Dmax;
kx = v';
ky = ones(nn,nn).*v;

kz =  repmat(ky(:),1,nn);
ky =  repmat(ky,nn,1);
kx =  repmat(kx,nn,1);
%
mref =(n+1)*(ky-1)+kx +(n+1)*(n+1)*(kz-1);
mref1=(n+1)*(ky-1)+kx +(n+1)*(n+1)*kz;
mref2=(n+1)*(ky-1)+kx +(n+1)*(n+1)*(kz+1);
%
macro_element=xn;

nvv(:,1) = mref(:);
nvv(:,2) = mref(:)+2;
nvv(:,3) = mref2(:)+2;
nvv(:,4) = mref2(:) ;

nvv(:,6) = mref(:)+2*n+4;
nvv(:,5) = mref(:)+2*n+2;

nvv(:,7) = mref2(:)+2*n+4;
nvv(:,8) = mref2(:)+2*n+2;

% front face
nvv(:,9) =  mref(:)+1;
nvv(:,10)=  mref(:)+n+3;
nvv(:,11) = mref(:)+2*n+3;
nvv(:,12)=  mref(:)+n+1;
nvv(:,13)=  mref(:)+n+2;

% central face
nvv(:,14) = mref1(:)+1;
nvv(:,15) = mref1(:)+2;
nvv(:,16)=  mref1(:)+n+3;
nvv(:,17) = mref1(:)+2*n+4;
nvv(:,18) = mref1(:)+2*n+3;
nvv(:,19) = mref1(:)+2*n+2;
nvv(:,20)=  mref1(:)+n+1;
nvv(:,21) = mref1(:);
nvv(:,22)=  mref1(:)+n+2;

% back face
nvv(:,23) = mref2(:)+1;
nvv(:,24)=  mref2(:)+n+3;
nvv(:,25) = mref2(:)+2*n+3;
nvv(:,26)=  mref2(:)+n+1;
nvv(:,27)=  mref2(:)+n+2;

mv(macro_element,1:27)=nvv(:,1:27);
%
macro_element = size(mv,1);

%
% down face: z=-Dmax
k1=find( xyz(:,3)==left );
e1=[];
% right face: x=Dmax
k2=find( xyz(:,1)==Dmax & xyz(:,2)<Dmax & xyz(:,2) >left & xyz(:,3)>left & xyz(:,3) <Dmax);
e2=[];
% top face: z=Dmax
k3=find( xyz(:,3)==Dmax);% & xyz(:,2)>-Dmax & xyz(:,2)<=Dmax & xyz(:,1)>-Dmax & xyz(:,1)<Dmax);
e3=[];
% left face: x=-Dmax
k4=find( xyz(:,1)==left & xyz(:,2)<Dmax  & xyz(:,2) >left  & xyz(:,3)<Dmax   & xyz(:,3) >left );
e4=[];
% front face: y=-Dmax
k5=find( xyz(:,2)==left & xyz(:,3)>left & xyz(:,3)<Dmax & xyz(:,1)<=Dmax  & xyz(:,1) >=left );
e5=[];
% back face: y=Dmax
k6=find( xyz(:,2)==Dmax & xyz(:,3)<Dmax   & xyz(:,3) >left & xyz(:,1)<=Dmax  & xyz(:,1) >=left);
e6=[];

% the following code is needed for plotting
for k=1:macro_element
    if any(mv(k,14)==k1)
        e1=[e1,k];
    end
    
    if any(mv(k,16)==k2)
        e2=[e2,k];
    end
    
    if any(mv(k,18)==k3)
        e3=[e3,k];
    end
    
    if any(mv(k,20)==k4)
        e4=[e4,k];
    end
    
    if any(mv(k,13)==k5)
        e5=[e5,k];
    end
    
    if any(mv(k,27)==k6)
        e6=[e6,k];
    end
    
end
ef1=ones(size(e1));
ef2=2*ones(size(e2));
ef3=3*ones(size(e3));
ef4=4*ones(size(e4));
ef5=5*ones(size(e5));
ef6=6*ones(size(e6));
mbound3D=[e1',ef1';e2',ef2';e3',ef3';e4',ef4';e5',ef5';e6',ef6'];

bound3D=sort([k1;k2;k3;k4;k5;k6]);

gohome
cd datafiles
save cube_grid.mat mv xyz bound3D mbound3D grid_type x y z n nc

return
