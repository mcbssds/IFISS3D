function bore_domain(Dmax,nc)
%BORE_DOMAIN borehole domain Q2-element grid generator
%   bore_domain(Dmax,nc);
%   inputs: 
%   Dmax   limits of domain [-Dmax, Dmax] 
%   nc     grid parameter 15*nc x 15*nc x 15*nc stretched grid
%
% grid defining data is saved to the file: bore_grid.mat
% IFISS function: DJS; 29 July 2022
% Copyright (c) G. Papanikos, C.E. Powell, D.J. Silvester

bd = 1e-2; % borehole semi-width
if nc<1,
    error('illegal parameter choice, try again.'),
end
if Dmax<0.1 || Dmax ==0
    error('domain is too small, try again.')
end

left =-Dmax;
nnx  = 5*nc;

% compute (x,y) coordinates of vertices
x_pos = subint_expand(Dmax,Dmax/2,bd,2*nnx,nnx); 
x_neg = -x_pos;
x = [x_neg(end:-1:1);0;x_pos];

% plot to illustrate borehole resolution
figure
plot(x,0*x,'bo'), axis('equal')
title('Location of mesh points in each direction')

nx = 6*nnx+3;
nz = nx;
z = x;
y = x;
n = nx-1;
np = n/2;
ny = length(y);
nvtx = nx*ny*nz;

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

% down face  z = -Dmax
k1=find( abs(xyz(:,3)-left)<10*eps );
e1=[];
% right face x = Dmax
k2=find( abs(xyz(:,1)-Dmax)<10*eps & xyz(:,2)<Dmax & xyz(:,2) >left & xyz(:,3)>left & xyz(:,3) <Dmax);
e2=[];
% top face  z = Dmax
k3=find( abs(xyz(:,3)-Dmax)<10*eps);% & xyz(:,2)>-Dmax & xyz(:,2)<=Dmax & xyz(:,1)>-Dmax & xyz(:,1)<Dmax);
e3=[];
% left face   x = -Dmax 
k4=find( abs(xyz(:,1)-left)<10*eps & xyz(:,2)<Dmax  & xyz(:,2) >left  & xyz(:,3)<Dmax   & xyz(:,3) >left );
e4=[];
% front face  y = -Dmax
k5=find( abs(xyz(:,2)-left)<10*eps & xyz(:,3)>left & xyz(:,3)<Dmax & xyz(:,1)<=Dmax  & xyz(:,1) >=left );
e5=[];
% back face   y = Dmax
k6=find( abs(xyz(:,2)-Dmax)<10*eps & xyz(:,3)<Dmax   & xyz(:,3) >left & xyz(:,1)<=Dmax  & xyz(:,1) >=left);
e6=[];

% Borehole values 
k7=find( abs(xyz(:,3)-bd)<10*eps & xyz(:,1)>-bd &  xyz(:,1)<bd & xyz(:,2)<Dmax &  xyz(:,2)>=0);
e7=[];
k8=find( abs(xyz(:,1)-bd)<10*eps & xyz(:,2)<Dmax &  xyz(:,2)>=0 & xyz(:,3)>-bd & xyz(:,3)< bd);
e8=[];
k9=find( abs(xyz(:,3)+bd)<10*eps & xyz(:,2)<Dmax &  xyz(:,2)>=0 & xyz(:,1)>-bd & xyz(:,1)< bd);
e9=[];
k10=find( abs(xyz(:,1)+bd)<10*eps & xyz(:,2)<Dmax &  xyz(:,2)>=0 & xyz(:,3)>-bd & xyz(:,3)< bd);
e10=[];
k11 = find (xyz(:,2) == 0 & xyz(:,3)>-bd & xyz(:,3)< bd & xyz(:,1)>-bd & xyz(:,1)< bd);

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
ef1=  ones(size(e1));
ef2=2*ones(size(e2));
ef3=3*ones(size(e3));
ef4=4*ones(size(e4));
ef5=5*ones(size(e5));
ef6=6*ones(size(e6));
mbound3D=[e1',ef1';e2',ef2';e3',ef3';e4',ef4';e5',ef5';e6',ef6'];

bound3D=sort([k1;k2;k3;k4;k5;k6;k7;k8;k9;k10;k11]);

gohome
cd datafiles
save -v7.3 bore_grid.mat mv xyz bound3D mbound3D x y z n nc

return
