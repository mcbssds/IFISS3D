function ev=q1grid3D(xyz,mv,bound3D);
%Q1GRID3D trilinear element grid generator
%   ev=q1grid3D(xyz,mv,bound3D);
%   inputs:
%          xyz        vertex coordinates
%          mv         Q2 macroelement mapping matrix
%          bound3D    boundary vertex vector
%   output:
%          ev         element mapping matrix
% IFISS function: GP; 9 June 2022.
% Copyright (c)   G.Papanikos, C.E. Powell, D.J. Silvester

xx=xyz(:,1); yy=xyz(:,2); zz=xyz(:,3); nvtx=length(xx);

mel = length(mv(:,1)); nel = 8*mel;
% ev=zeros(nel,8);
%
%% create elements from the macroelements 
k =1:mel;
ke=8*k-7;
ev(ke,1)=mv(k,1);
ev(ke,2)=mv(k,9);
ev(ke,3)=mv(k,14);
ev(ke,4)=mv(k,21);
ev(ke,5)=mv(k,12);
ev(ke,6)=mv(k,13);
ev(ke,7)=mv(k,22);
ev(ke,8)=mv(k,20);

% second element
ke=8*k-6;
ev(ke,1)=mv(k,9);
ev(ke,2)=mv(k,2);
ev(ke,3)=mv(k,15);
ev(ke,4)=mv(k,14);
ev(ke,5)=mv(k,13);
ev(ke,6)=mv(k,10);
ev(ke,7)=mv(k,16);
ev(ke,8)=mv(k,22);

%third element 
ke=8*k-5;
ev(ke,1)=mv(k,14);
ev(ke,2)=mv(k,15);
ev(ke,3)=mv(k,3);
ev(ke,4)=mv(k,23);
ev(ke,5)=mv(k,22);
ev(ke,6)=mv(k,16);
ev(ke,7)=mv(k,24);
ev(ke,8)=mv(k,27);

%fourth element 
ke=8*k-4;
ev(ke,1)=mv(k,21);
ev(ke,2)=mv(k,14);
ev(ke,3)=mv(k,23);
ev(ke,4)=mv(k,4);
ev(ke,5)=mv(k,20);
ev(ke,6)=mv(k,22);
ev(ke,7)=mv(k,27);
ev(ke,8)=mv(k,26);

%fifth element
ke=8*k-3;
ev(ke,1)=mv(k,12);
ev(ke,2)=mv(k,13);
ev(ke,3)=mv(k,22);
ev(ke,4)=mv(k,20);
ev(ke,5)=mv(k,5);
ev(ke,6)=mv(k,11);
ev(ke,7)=mv(k,18);
ev(ke,8)=mv(k,19);

% sixth element
ke=8*k-2;
ev(ke,1)=mv(k,13);
ev(ke,2)=mv(k,10);
ev(ke,3)=mv(k,16);
ev(ke,4)=mv(k,22);
ev(ke,5)=mv(k,11);
ev(ke,6)=mv(k,6);
ev(ke,7)=mv(k,17);
ev(ke,8)=mv(k,18);

%seventh element 
ke=8*k-1;
ev(ke,1)=mv(k,22);
ev(ke,2)=mv(k,16);
ev(ke,3)=mv(k,24);
ev(ke,4)=mv(k,27);
ev(ke,5)=mv(k,18);
ev(ke,6)=mv(k,17);
ev(ke,7)=mv(k,7);
ev(ke,8)=mv(k,25);

%eighth element 
ke=8*k;
ev(ke,1)=mv(k,20);
ev(ke,2)=mv(k,22);
ev(ke,3)=mv(k,27);
ev(ke,4)=mv(k,26);
ev(ke,5)=mv(k,19);
ev(ke,6)=mv(k,18);
ev(ke,7)=mv(k,25);
ev(ke,8)=mv(k,8);

return
