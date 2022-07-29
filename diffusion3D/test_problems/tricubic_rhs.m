function f = specific_rhs3D(x,y,z,nel,ngpt)
% RHS forcing function
%   f = specific_rhs3D(x,y,z,nel)
%   input
%          x          x coordinate vector
%          y          y coordinate vector 
%          z          y coordinate vector 
%          nel        number of elements  
% IFISS function: DJS; 27 July 2022.
% Copyright (c)  2022  G.Papanikos,  C.E. Powell, D.J. Silvester

 f = 6*x.*y.*z.*(y.^2 - 1).*(z.^2 - 1) + ...
     6*x.*y.*z.*(x.^2 - 1).*(z.^2 - 1) + ...
     6*x.*y.*z.*(x.^2 - 1).*(y.^2 - 1);     % tricubic solution
return
