function [psi,dpsids,dpsidt,dpsidl] = qshape3D(s,t,l)
%QSHAPE3D evaluates triquadratic shape functions 
%   [psi,dpsids,dpsidt] = qshape3D(s,t,l);
%   input
%          s         x coordinate   
%          t         y coordinate
%          l         z coordinate 
%   output
%          psi        shape function
%          dpsids     x derivative of psi
%          dpsidt     y derivative of psi
%          dpsidl     z derivative of psi
% IFISS function: GP; 9 June 2022.
% Copyright (c)  2022  G.Papanikos,  C.E. Powell, D.J. Silvester
%
% --------------------- Node numbering -----------------------------------
%
%
%

% one dimensional shape functions      
      ellx(1) = 0.5*s*(s-1);  elly(1) = 0.5*t*(t-1); ellz(1) = 0.5*l*(l-1);
      ellx(2) = 1-(s*s);      elly(2) = 1-(t*t);     ellz(2) = 1-(l*l);
	  ellx(3) = 0.5*s*(s+1);  elly(3) = 0.5*t*(t+1); ellz(3) = 0.5*l*(l+1);
      
      dellx(1) = s-0.5;       delly(1) = t-0.5;      dellz(1) = l-0.5;
	  dellx(2) = -2*s;        delly(2) = -2*t;       dellz(2) = -2*l;
	  dellx(3) = s+0.5;       delly(3) = t+0.5;      dellz(3) = l+0.5;
      
    
%3D dimensional shape functions--- derivative w.t.r s ------------------- derivative w.r.t t -------------------  derivative w.r.t l
 psi(1) = ellx(1)*elly(1)*ellz(1);  dpsids(1) = dellx(1)*elly(1)*ellz(1);  dpsidt(1) = ellx(1)*delly(1)*ellz(1);   dpsidl(1) = ellx(1)*elly(1)*dellz(1);  %point (-1,-1,-1) 1
 psi(2) = ellx(3)*elly(1)*ellz(1);  dpsids(2) = dellx(3)*elly(1)*ellz(1);  dpsidt(2) = ellx(3)*delly(1)*ellz(1);   dpsidl(2) = ellx(3)*elly(1)*dellz(1);  %point ( 1,-1,-1) 2
 psi(3) = ellx(3)*elly(3)*ellz(1);  dpsids(3) = dellx(3)*elly(3)*ellz(1);  dpsidt(3) = ellx(3)*delly(3)*ellz(1);   dpsidl(3) = ellx(3)*elly(3)*dellz(1);  %point ( 1, 1,-1) 3  
 psi(4) = ellx(1)*elly(3)*ellz(1);  dpsids(4) = dellx(1)*elly(3)*ellz(1);  dpsidt(4) = ellx(1)*delly(3)*ellz(1);   dpsidl(4) = ellx(1)*elly(3)*dellz(1);  %point (-1, 1,-1) 4
 psi(5) = ellx(1)*elly(1)*ellz(3);  dpsids(5) = dellx(1)*elly(1)*ellz(3);  dpsidt(5) = ellx(1)*delly(1)*ellz(3);   dpsidl(5) = ellx(1)*elly(1)*dellz(3);  %point (-1,-1, 1) 5
 psi(6) = ellx(3)*elly(1)*ellz(3);  dpsids(6) = dellx(3)*elly(1)*ellz(3);  dpsidt(6) = ellx(3)*delly(1)*ellz(3);   dpsidl(6) = ellx(3)*elly(1)*dellz(3);  %point ( 1,-1, 1) 6
 psi(7) = ellx(3)*elly(3)*ellz(3);  dpsids(7) = dellx(3)*elly(3)*ellz(3);  dpsidt(7) = ellx(3)*delly(3)*ellz(3);   dpsidl(7) = ellx(3)*elly(3)*dellz(3);  %point ( 1, 1, 1) 7  
 psi(8) = ellx(1)*elly(3)*ellz(3);  dpsids(8) = dellx(1)*elly(3)*ellz(3);  dpsidt(8) = ellx(1)*delly(3)*ellz(3);   dpsidl(8) = ellx(1)*elly(3)*dellz(3);  %point (-1, 1, 1) 8
 

 psi(9)  = ellx(2)*elly(1)*ellz(1); dpsids(9)  = dellx(2)*elly(1)*ellz(1); dpsidt(9)  = ellx(2)*delly(1)*ellz(1);  dpsidl(9)  = ellx(2)*elly(1)*dellz(1); %point ( 0,-1,-1) 9  
 psi(10) = ellx(3)*elly(1)*ellz(2); dpsids(10) = dellx(3)*elly(1)*ellz(2); dpsidt(10) = ellx(3)*delly(1)*ellz(2);  dpsidl(10) = ellx(3)*elly(1)*dellz(2); %point ( 1,-1, 0) 10  
 psi(11) = ellx(2)*elly(1)*ellz(3); dpsids(11) = dellx(2)*elly(1)*ellz(3); dpsidt(11) = ellx(2)*delly(1)*ellz(3);  dpsidl(11) = ellx(2)*elly(1)*dellz(3); %point ( 0,-1, 1) 11 
 psi(12) = ellx(1)*elly(1)*ellz(2); dpsids(12) = dellx(1)*elly(1)*ellz(2); dpsidt(12) = ellx(1)*delly(1)*ellz(2);  dpsidl(12) = ellx(1)*elly(1)*dellz(2); %point (-1,-1, 0) 12  
 psi(13) = ellx(2)*elly(1)*ellz(2); dpsids(13) = dellx(2)*elly(1)*ellz(2); dpsidt(13) = ellx(2)*delly(1)*ellz(2);  dpsidl(13) = ellx(2)*elly(1)*dellz(2); %point ( 0,-1, 0) 13 

 psi(14) = ellx(2)*elly(2)*ellz(1); dpsids(14) = dellx(2)*elly(2)*ellz(1); dpsidt(14) = ellx(2)*delly(2)*ellz(1);  dpsidl(14) = ellx(2)*elly(2)*dellz(1); %point ( 0, 0,-1) 14 
 psi(15) = ellx(3)*elly(2)*ellz(1); dpsids(15) = dellx(3)*elly(2)*ellz(1); dpsidt(15) = ellx(3)*delly(2)*ellz(1);  dpsidl(15) = ellx(3)*elly(2)*dellz(1); %point ( 1, 0,-1) 15  
 psi(16) = ellx(3)*elly(2)*ellz(2); dpsids(16) = dellx(3)*elly(2)*ellz(2); dpsidt(16) = ellx(3)*delly(2)*ellz(2);  dpsidl(16) = ellx(3)*elly(2)*dellz(2); %point ( 1, 0, 0) 16  
 psi(17) = ellx(3)*elly(2)*ellz(3); dpsids(17) = dellx(3)*elly(2)*ellz(3); dpsidt(17) = ellx(3)*delly(2)*ellz(3);  dpsidl(17) = ellx(3)*elly(2)*dellz(3); %point ( 1, 0, 1) 17   
 psi(18) = ellx(2)*elly(2)*ellz(3); dpsids(18) = dellx(2)*elly(2)*ellz(3); dpsidt(18) = ellx(2)*delly(2)*ellz(3);  dpsidl(18) = ellx(2)*elly(2)*dellz(3); %point ( 0, 0, 1) 18 
 psi(19) = ellx(1)*elly(2)*ellz(3); dpsids(19) = dellx(1)*elly(2)*ellz(3); dpsidt(19) = ellx(1)*delly(2)*ellz(3);  dpsidl(19) = ellx(1)*elly(2)*dellz(3); %point (-1, 0, 1) 19  
 psi(20) = ellx(1)*elly(2)*ellz(2); dpsids(20) = dellx(1)*elly(2)*ellz(2); dpsidt(20) = ellx(1)*delly(2)*ellz(2);  dpsidl(20) = ellx(1)*elly(2)*dellz(2); %point (-1, 0 ,0) 20  
 psi(21) = ellx(1)*elly(2)*ellz(1); dpsids(21) = dellx(1)*elly(2)*ellz(1); dpsidt(21) = ellx(1)*delly(2)*ellz(1);  dpsidl(21) = ellx(1)*elly(2)*dellz(1); %point (-1, 0,-1) 21
 psi(22) = ellx(2)*elly(2)*ellz(2); dpsids(22) = dellx(2)*elly(2)*ellz(2); dpsidt(22) = ellx(2)*delly(2)*ellz(2);  dpsidl(22) = ellx(2)*elly(2)*dellz(2); %point ( 0, 0, 0) 22 

 
 psi(23) = ellx(2)*elly(3)*ellz(1); dpsids(23) = dellx(2)*elly(3)*ellz(1); dpsidt(23) = ellx(2)*delly(3)*ellz(1);  dpsidl(23) = ellx(2)*elly(3)*dellz(1); %point ( 0, 1,-1) 23  
 psi(24) = ellx(3)*elly(3)*ellz(2); dpsids(24) = dellx(3)*elly(3)*ellz(2); dpsidt(24) = ellx(3)*delly(3)*ellz(2);  dpsidl(24) = ellx(3)*elly(3)*dellz(2); %point ( 1, 1, 0) 24  
 psi(25) = ellx(2)*elly(3)*ellz(3); dpsids(25) = dellx(2)*elly(3)*ellz(3); dpsidt(25) = ellx(2)*delly(3)*ellz(3);  dpsidl(25) = ellx(2)*elly(3)*dellz(3); %point ( 0, 1, 1) 25  
 psi(26) = ellx(1)*elly(3)*ellz(2); dpsids(26) = dellx(1)*elly(3)*ellz(2); dpsidt(26) = ellx(1)*delly(3)*ellz(2);  dpsidl(26) = ellx(1)*elly(3)*dellz(2); %point (-1, 1, 0) 26  
 psi(27) = ellx(2)*elly(3)*ellz(2); dpsids(27) = dellx(2)*elly(3)*ellz(2); dpsidt(27) = ellx(2)*delly(3)*ellz(2);  dpsidl(27) = ellx(2)*elly(3)*dellz(2); %point ( 0, 1, 0) 27  



 
 





 end