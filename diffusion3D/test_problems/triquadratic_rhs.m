function f = specific_rhs3D(x,y,z,nel)
%TRIQUADRATIC_RHS generates specific forcing function f
%   f = specific_rhs3D(x,y,z,nel)
%   inputs:
%          x          x coordinate vector
%          y          y coordinate vector 
%          z          y coordinate vector 
%          nel        number of elements  
% IFISS function: DJS; 4 September 2022.
% Copyright (c)  2022  G.Papanikos,  C.E. Powell, D.J. Silvester

% specific f corresponding to triquadratic solution
 f= 2*(x.^2 - 1).*(y.^2 - 1) + 2*(x.^2 - 1).*(z.^2 - 1) + ...
    2*(y.^2 -1).*(z.^2 - 1);         
return
