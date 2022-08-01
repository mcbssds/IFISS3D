 function [phi,dphids,dphidt,dphidl] = shape3D(s,t,l)
%SHAPE3D evaluates trilinear shape functions 
%   [phi,dphids,dphidt,dphidl] = shape3D(s,t,l);
%   inputs:
%          s         x coordinate   
%          t         y coordinate
%          l         z coordinate 
%   outputs:
%          phi        shape function
%          dphids     x derivative of phi
%          dphidt     y derivative of phi
%          dphidl     z derivative of phi 
%
%                                         l
%           8 _ _ _ _7                    |                      
%          /|       /|                    |
%         / |      / |                    |_ _ _ _s          
%        5_ |_  _6_ _|                   /                   
%        |  /4    |  /3                 /
%        | /      | /                  /t
%        |/_ _ _ _|/  
%        1       2
%
% IFISS function: GP; 9 June 2022.
% Copyright (c)  2022  G.Papanikos,  C.E. Powell, D.J. Silvester
      one = 1.0e0;
%     trilinear basis functions 
      phi(1,:) = 1/8 *(one-s).*(one-t).*(one-l);  %(-1,-1,-1) 1
      phi(2,:) = 1/8 *(one+s).*(one-t).*(one-l);  %( 1,-1,-1) 2
      phi(3,:) = 1/8 *(one+s).*(one+t).*(one-l);  %( 1, 1,-1) 3
      phi(4,:) = 1/8 *(one-s).*(one+t).*(one-l);  %(-1, 1,-1) 4
      phi(5,:) = 1/8 *(one-s).*(one-t).*(one+l);  %(-1,-1, 1) 5
      phi(6,:) = 1/8 *(one+s).*(one-t).*(one+l);  %( 1,-1, 1) 6
      phi(7,:) = 1/8 *(one+s).*(one+t).*(one+l);  %( 1, 1, 1) 7
      phi(8,:) = 1/8 *(one-s).*(one+t).*(one+l);  %(-1, 1, 1) 8
 
      % derivatives with respect to s
      dphids(1,:) = -1/8 *(one-t).*(one-l);
      dphids(2,:) =  1/8 *(one-t).*(one-l);
      dphids(3,:) =  1/8 *(one+t).*(one-l);
      dphids(4,:) = -1/8 *(one+t).*(one-l);
      dphids(5,:) = -1/8 *(one-t).*(one+l);
      dphids(6,:) =  1/8 *(one-t).*(one+l);
      dphids(7,:) =  1/8 *(one+t).*(one+l);
      dphids(8,:) = -1/8 *(one+t).*(one+l);
      
      % derivatives with respect to t
      dphidt(1,:) = -1/8 *(one-s).*(one-l);
      dphidt(2,:) = -1/8 *(one+s).*(one-l);
      dphidt(3,:) =  1/8 *(one+s).*(one-l);
      dphidt(4,:) =  1/8 *(one-s).*(one-l);
      dphidt(5,:) = -1/8 *(one-s).*(one+l);
      dphidt(6,:) = -1/8 *(one+s).*(one+l);
      dphidt(7,:) =  1/8 *(one+s).*(one+l);
      dphidt(8,:) =  1/8 *(one-s).*(one+l);
      
      % derivatives with respect to l
      dphidl(1,:) = -1/8 *(one-s).*(one-t);
      dphidl(2,:) = -1/8 *(one+s).*(one-t);
      dphidl(3,:) = -1/8 *(one+s).*(one+t);
      dphidl(4,:) = -1/8 *(one-s).*(one+t);
      dphidl(5,:) =  1/8 *(one-s).*(one-t);
      dphidl(6,:) =  1/8 *(one+s).*(one-t);
      dphidl(7,:) =  1/8 *(one+s).*(one+t);
      dphidl(8,:) =  1/8 *(one-s).*(one+t);
      
      return
