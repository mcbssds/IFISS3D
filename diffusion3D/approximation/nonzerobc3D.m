function [a,f] = nonzerobc3D(a,f,xyz,bound3D)
%NONZEROBC3D imposes Dirichlet boundary condition 
%   [a,f] = nonzerobc3D(A,f,xyz,bound3D);
%   input
%          a          stiffness matrix
%          f          rhs vector
%          xyz         vertex coordinate vector  
%          bound3D      boundary vertex vector
%   output
%          Agal       stiffness matrix
%          fgal       rhs vector
%
%   calls function specific_bc
% IFISS function: GP; 9 June 2022; .
% Copyright (c) 2022 G. Papanikos, C.E. Powell, D.J. Silvester

nvtx = length(f); nbd=length(bound3D);
null_col=sparse(nvtx,nbd); null_row=sparse(nbd,nvtx);
Ax=a(1:nvtx,1:nvtx);
fx=f(1:nvtx);
%% set boundary condition
xbd=xyz(bound3D,1); ybd=xyz(bound3D,2); zbd=xyz(bound3D,3);
bc=specific_bc3D(xbd,ybd,zbd); 
fx = fx - Ax(:,bound3D)*bc;
dA=zeros(nvtx,1); dA(bound3D)=ones(nbd,1);
Ax(:,bound3D)=null_col;  Ax(bound3D,:)=null_row;   
Ax=Ax+spdiags(dA,0,nvtx,nvtx);  fx(bound3D)=bc; 
a=Ax; f=fx;
return
