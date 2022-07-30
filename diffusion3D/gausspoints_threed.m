function [s,t,l,w]=gausspoints_threed(oneg,onew) 
%GAUSSPOINTS_THREED constructs 3D tensor product Gaussian quadrature rule
% [s,t,l,w] = gausspoints_threed(oneg,onew);
%  inputs:
%       oneg    Gauss points for 1D rule in (-1,1)
%       onew    weights of 1D rule
%  outputs:
%       s       x-coordinates of the 3D Gauss points in (-1,1)^{3}
%       t       y-coordinates of the 3D Gauss points in (-1,1)^{3}
%       l       z-coordinates of the 3D Gauss points in (-1,1)^{3}
%       w       weights of the tensor product 3D rule
% IFISS function: CP; 30th July 2022.
% Copyright (c)  2022  G.Papanikos, C.E. Powell, D.J. Silvester

ng=max(size(oneg)); 
s=zeros(ng*ng*ng,1);
t=zeros(ng*ng*ng,1);
l=zeros(ng*ng*ng,1);
w=zeros(ng*ng*ng,1);
for i=1:ng
    for j=1:ng
        for k =1:ng
        s(i+(j-1)*ng+ng*ng*(k-1))=oneg(k);
        t(i+(j-1)*ng+ng*ng*(k-1))=oneg(j);
        l(i+(j-1)*ng+ng*ng*(k-1))=oneg(i);
        w(i+(j-1)*ng+ng*ng*(k-1))=onew(i)*onew(j)*onew(k);
        end
    end
end
return
