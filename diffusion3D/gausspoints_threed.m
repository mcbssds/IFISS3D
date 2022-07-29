function [s,t,l,w]=gausspoints_threed(oneg,onew) 
%GAUSSPOINTS_THREED constructs tensor product Gauss Point Rule
% [s,t,l,w] = gausspoints_threed(oneg,onew);
%  input
%   oneg    1D Gaussian Points in (-1,1)
%   onew    Weights in 1D
%  output
%   s       x-coordinate of the Gaussian Points
%   t       y-coordinate of the Gaussian Points
%   l       z-coordinate of the Gaussian Points
%   w       weights in 3D
% IFISS scriptfile: GP; 9 June 2022.
% Copyright (c)  2022  G.Papanikos,  C.E. Powell, D.J. Silvester

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
