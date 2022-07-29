function f = specific_rhs3D(x,y,z,nel,ngpt)
%unit_rhs   unit RHS forcing function
%   f = specific_rhs(x,y,z,nel)
%   input
%          x          x coordinate vector
%          y          y coordinate vector 
%          z          y coordinate vector 
%          nel        number of elements  
% IFISS function: GP; 9 June 2022.
% Copyright (c)  2022  G.Papanikos,  C.E. Powell, D.J. Silvester

 f=ones(nel,ngpt);
 
return