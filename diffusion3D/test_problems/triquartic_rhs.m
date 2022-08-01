function f = specific_rhs3D(x,y,z,nel,ngpt)
%TRIQUARTIC_RHS evaluates specific forcing function f
%   f = specific_rhs3D(x,y,z,nel)
%   inputs:
%          x          x coordinate vector
%          y          y coordinate vector 
%          z          y coordinate vector 
%          nel        number of elements  
% IFISS function: GP; 9 June 2022.
% Copyright (c)  2022  G.Papanikos,  C.E. Powell, D.J. Silvester

% specific f corresponding to triquartic solution
f = 12*x.^2.*(y.^4 - 1).*(z.^4 - 1) + 12*y.^2.*(x.^4 - 1).*(z.^4 - 1) + 12*z.^2.*(x.^4 - 1).*(y.^4 - 1); 
 
return
