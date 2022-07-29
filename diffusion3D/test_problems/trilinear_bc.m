function bc = specific_bc3D(xbd,ybd,zbd)
%bilinear_bc 
%   bc = specific_bc3D(xbd,ybd,zbd);
%   input
%          xbd          x boundary coordinate vector
%          ybd          y boundary coordinate vector 
%          zbd          z boundary coordinate vector 
% IFISS function: GP; 9 June 2022.
% Copyright (c)  2022  G.Papanikos,  C.E. Powell, D.J. Silvester


 indx  =find(xbd== 1);
 indx2  =find(xbd== -1);
% 
 indy  =find(ybd == 1);
 indy2  =find(ybd == -1);
% 
 indz  =find(zbd == 1);
 indz2  =find(zbd == -1);
 
 bc(indx) =0;
 bc(indx2) = 2.*(ybd(indx2)-1).*(zbd(indx2)-1);
% 
% 
 bc(indy) =0;
 bc(indy2) = 2.*(xbd(indy2)-1).*(zbd(indy2)-1);
 
 bc(indz) =0;
 bc(indz2) = 2.*(xbd(indz2)-1).*(ybd(indz2)-1);
bc =bc';

return