function stair_domain
%STAIR_DOMAIN stair-shaped domain Q2 grid generator
%   stair_domain;
% 
% grid defining data is saved to the file: stair_grid.mat
% IFISS function: DJS; 29 July 2022.
% Copyright (c) 2022 G. Papanikos,  C.E. Powell, D.J. Silvester

Dmax = 1;
fprintf('\n\nGrid generation for a stair-shaped domain.\n')
nc=default('grid parameter: 3 for underlying 8x8x8 grid (default is 4 for 16x16x16 grid)',4);
if nc<2, error('illegal parameter choice, try again.'), end
n=2^nc; np=n/(Dmax*2); nq=n/(Dmax*4);

% compute (x,y,z) coordinates of vertices
% y-direction
yy=[1/np:1/np:Dmax];
ypos=[0,yy];
yneg=-yy(length(yy):-1:1);
y=[yneg,ypos];
% x-direction
xx=yy;
xpos=[0,xx];
xneg=-xx(length(xx):-1:1);
x=[xneg,xpos];
%
zz=yy;
zpos=[0,zz];
zneg=-zz(length(zz):-1:1);
z=[zneg,zpos];

% compute triquadratic element coordinates
nvtx=(Dmax*nq+1)*(Dmax*nq)*(Dmax*n);
% negative x-values
[Xpos,Yneg,Z]=meshgrid(xpos,yneg,z);
xx=reshape(Yneg,Dmax*np*(Dmax*np+1)*(2*Dmax*np+1),1);
yy=reshape(Z,Dmax*np*(Dmax*np+1)*(2*Dmax*np+1),1);
zz=reshape(Xpos,Dmax*np*(Dmax*np+1)*(2*Dmax*np+1),1);

xyzleft=[xx(:),yy(:),zz(:)];
kz =1;
v  = 1:2:nq*Dmax+nq;
nn = length(v);
xn = 1:nq*nq*Dmax;
kx = v;
ky = ones(nn,nn).*v;

for k=1:np*Dmax
    mref=  Dmax*np*(ky-1)+kx' +(np+1)*np*(kz-1);
    mref1= Dmax*np*(ky-1)+kx' +(np+1)*np*kz;
    mref2= Dmax*np*(ky-1)+kx' +(np+1)*np*(kz+1);

    macro_element=xn;

    nvv(:,1) = mref(:);
    nvv(:,2) = mref(:)+2;
    nvv(:,3) = mref2(:)+2;
    nvv(:,4) = mref2(:) ;

    nvv(:,6) = mref(:)+2*Dmax*np+2;
    nvv(:,5) = mref(:)+2*Dmax*np;

    nvv(:,7) = mref2(:)+2*Dmax*np+2;
    nvv(:,8) = mref2(:)+2*Dmax*np;

    %face front

    nvv(:,9) =  mref(:)+1;
    nvv(:,10)=  mref(:)+Dmax*np+2;
    nvv(:,11) = mref(:)+2*Dmax*np+1;
    nvv(:,12)=  mref(:)+Dmax*np;
    nvv(:,13)=  mref(:)+Dmax*np+1;


    % central face

    nvv(:,14) = mref1(:)+1;
    nvv(:,15) = mref1(:)+2;
    nvv(:,16)=  mref1(:)+Dmax*np+2;
    nvv(:,17) = mref1(:)+2*Dmax*np+2;
    nvv(:,18) = mref1(:)+2*Dmax*np+1;
    nvv(:,19) = mref1(:)+2*Dmax*np;
    nvv(:,20)=  mref1(:)+Dmax*np;
    nvv(:,21) = mref1(:);
    nvv(:,22)=  mref1(:)+Dmax*np+1;


    % back face
    nvv(:,23) = mref2(:)+1;
    nvv(:,24)=  mref2(:)+Dmax*np+2;
    nvv(:,25) = mref2(:)+2*Dmax*np+1;
    nvv(:,26)=  mref2(:)+Dmax*np;
    nvv(:,27)=  mref2(:)+Dmax*np+1;


    mv(macro_element,1:27)=nvv(:,1:27);
    xn = xn + nq*nq;
    kz =kz+2;
end
   
[X,Ypos,Z]=meshgrid(x,ypos,z);
xx=reshape(Ypos,(Dmax*np+1)*(n+1)*(n+1),1);
yy=reshape(Z,(Dmax*np+1)*(n+1)*(n+1),1);
zz =reshape(X,(Dmax*np+1)*(n+1)*(n+1),1);
xyzright=[xx(:),yy(:),zz(:)];
xyz=[xyzleft;xyzright]; 

% correction along the internal boundary
ky=1;
kx= 1;
kz =1;
l =0;
for k= 1:Dmax*np
mref = (Dmax*np+1)*(ky-1)+kx +(np+1)*(n+1)*(kz-1)+length(xyzleft(:,1))+(np+1)*(np);
mref1= (Dmax*np+1)*(ky-1)+kx  +(np+1)*(n+1)*(kz)+length(xyzleft(:,1))+(np+1)*(np);
mref2= (Dmax*np+1)*(ky-1)+kx  +(np+1)*(n+1)*(kz+1)+length(xyzleft(:,1))+(np+1)*(np);
      
for mel=Dmax*nq*(l+1):Dmax*nq:Dmax*nq*Dmax*nq*k;
nvv=mv(mel,:);
nvv(2)  = mref;
nvv(6)  = mref+2*Dmax*np+2;
nvv(10) = mref+Dmax*np+1; 

nvv(15)  =  mref1;
nvv(17)  =  mref1+2*Dmax*np+2;
nvv(16)  =  mref1+Dmax*np+1; 

nvv(3)  =  mref2;
nvv(7)  =  mref2+2*Dmax*np+2;
nvv(24) =  mref2+Dmax*np+1; 

mv(mel,1:27)=nvv(1:27);
mref=mref+2*Dmax*np+2;
mref1=mref1+2*Dmax*np+2;
mref2=mref2+2*Dmax*np+2;
end
kz= kz+2;
l=l+nq;
end

clear nvv
ky = 1;
kz =1;
xn = 1:nq*Dmax;
xn = xn+length(mv(:,1));
kx = 1:2:nq*Dmax+nq;

for k=1:np*Dmax
for j=1:np*Dmax
      mref = (Dmax*np+1)*(ky-1)+kx +(np+1)*(n+1)*(kz-1)+length(xyzleft(:,1));
      mref1= (Dmax*np+1)*(ky-1)+kx  +(np+1)*(n+1)*(kz)+length(xyzleft(:,1));
      mref2= (Dmax*np+1)*(ky-1)+kx  +(np+1)*(n+1)*(kz+1)+length(xyzleft(:,1));
      
      macro_element=xn;
     
      nvv(:,1) = mref;
      nvv(:,2) = mref+2;
      nvv(:,3) = mref2+2;
      nvv(:,4) = mref2 ;
      
      nvv(:,6) = mref+2*Dmax*np+4;
      nvv(:,5) = mref+2*Dmax*np+2;
      
      nvv(:,7) = mref2+2*Dmax*np+4;
      nvv(:,8) = mref2+2*Dmax*np+2;
      
      % front face
      nvv(:,9) =  mref+1;
      nvv(:,10)=  mref+Dmax*np+3;
      nvv(:,11) = mref+2*Dmax*np+3;
      nvv(:,12)=  mref+Dmax*np+1;
      nvv(:,13)=  mref+Dmax*np+2;
           
      % central face
      nvv(:,14) = mref1+1;
      nvv(:,15) = mref1+2;
      nvv(:,16)=  mref1+Dmax*np+3;
      nvv(:,17) = mref1+2*Dmax*np+4;
      nvv(:,18) = mref1+2*Dmax*np+3;
      nvv(:,19) = mref1+2*Dmax*np+2;
      nvv(:,20)=  mref1+Dmax*np+1;
      nvv(:,21) = mref1;
      nvv(:,22)=  mref1+Dmax*np+2;
      
      % back face
      nvv(:,23) = mref2+1;
      nvv(:,24)=  mref2+Dmax*np+3;
      nvv(:,25) = mref2+2*Dmax*np+3;
      nvv(:,26)=  mref2+Dmax*np+1;
      nvv(:,27)=  mref2+Dmax*np+2;
           
      mv(macro_element,1:27)=nvv(:,1:27);
    
   ky = ky + 2; 
   xn = xn+nq;
end
   ky=1;
   kz =kz+2;
end
 macro_element = size(mv,1);
%
%% compute boundary vertices
% eight boundary edges 

k1 = find(xyz(:,3)==0 & xyz(:,1) <0 );
k2 =find( xyz(:,1)==0  & xyz(:,3)<=0 & xyz(:,3)>-Dmax);
% down face  z=-Dmax
k3 = find(xyz(:,3)==-Dmax );
% x=Dmax
k4=find( xyz(:,1)==Dmax & xyz(:,2)<Dmax & xyz(:,2) >-Dmax & xyz(:,3)>-Dmax & xyz(:,3) <Dmax);
% z=Dmax
k5=find( xyz(:,3)==Dmax & xyz(:,2)>-Dmax & xyz(:,2)<Dmax);
% x=-Dmax
k6=find( xyz(:,1)==-Dmax & xyz(:,3)>0 & xyz(:,3)<Dmax & xyz(:,2)>-Dmax & xyz(:,2)<Dmax);
% front face  y=-1
k71=find( xyz(:,2)==-Dmax & xyz(:,3)>0 );
k72= find(xyz(:,2)==-Dmax & xyz(:,3)<=0 & xyz(:,3)>-Dmax & xyz(:,1)>0);
k7 =[k71;k72];
% y =Dmax
k81=find( xyz(:,2)== Dmax & xyz(:,3)>0);
k82= find(xyz(:,2)== Dmax & xyz(:,3)<=0 & xyz(:,3)>-Dmax & xyz(:,1)>0);
k8 =[k81;k82];

e1=[];
e2=[];
e3=[];
e4=[];
e5=[];
e6=[];
e7=[];
e8=[];

for k=1:macro_element
        
    if any(mv(k,14)==k1), e1=[e1,k]; end,
    if any(mv(k,12)==k2), e2=[e2,k]; end,
    if any(mv(k,9)==k3),  e3=[e3,k]; end, 
    if any(mv(k,16)==k4), e4=[e4,k]; end,
    if any(mv(k,18)==k5), e5=[e5,k]; end, 
    if any(mv(k,20)==k6), e6=[e6,k]; end,
    if any(mv(k,13)==k7), e7=[e7,k]; end, 
    if any(mv(k,27)==k8), e8=[e8,k]; end,

end

ef1=ones(size(e1));
ef2=2*ones(size(e2));
ef3=3*ones(size(e3));
ef4=4*ones(size(e4));
ef5=5*ones(size(e5));
ef6=6*ones(size(e6));
ef7=7*ones(size(e7));
ef8=8*ones(size(e8));

bound3D=sort([k1;k2;k3;k4;k5;k6;k7;k8]);
mbound3D=[e1',ef1';e2',ef2';e3',ef3';e4',ef4';e5',ef5';e6',ef6';e7',ef7';e8',ef8'];
%mbound3D=[];

%
gohome 
cd datafiles
save stair_grid.mat mv xyz bound3D mbound3D x y z xyzleft xyzright nc
return
