function f = specific_rhs3D(x,y,z,nel)
%ZERO_RHS3D generates zero source function f=0
%   f = specific_rhs(x,y,z,nel)
%   inputs:
%          x          x coordinate vector
%          y          y coordinate vector
%          z          z coordinate vector
%          nel        number of elements
% IFISS function: DJS; 4 September 2022
% Copyright (c)  2022  G.Papanikos,  C.E. Powell, D.J. Silvester

f=0.*x;
return
