function bc = specific_bc3D(xbd,ybd,zbd)
%zero_bc   zero boundary condition
%   bc = specific_bc3D(xbd,ybd,zbd);
%   input
%          xbd          x boundary coordinate vector
%          ybd          y boundary coordinate vector
%          zbd          z boundary coordinate vector
% IFISS function: GP; 9 June 2022.
% Copyright (c)  2022  G.Papanikos,  C.E. Powell, D.J. Silvester
bc=zeros(size(xbd));
return