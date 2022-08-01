function f = specific_rhs3D(x,y,z,nel,ngpt)
%UNIT_RHS3D evaluates source function f=1
%   f = specific_rhs(x,y,z,nel)
%   inputs:
%          x          x coordinate vector
%          y          y coordinate vector 
%          z          y coordinate vector 
%          nel        number of elements  
%          ngpt       ???
% IFISS function: GP; 9 June 2022.
% Copyright (c)  2022  G.Papanikos,  C.E. Powell, D.J. Silvester

 f=ones(nel,ngpt);
 
return
