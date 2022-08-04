function [elerr_p] = diffpost_bc3D_Q1half(elerror,fez,xyz,ev,ebound3D,xl_m,yl_m,zl_m,xl_s,yl_s,zl_s)
%DIFFPOST_BC3D postprocesses local Poisson error estimator for full Q1(h/2)
%space
%   [elerr_p] = diffpost_bc3D_Q1half(elerror,fez,xyz,ev,ebound3D,xl_m,yl_m,zl_m,xl_s,yl_s,zl_s)
%   input
%          elerror       element error estimate (without BC imposition)
%          fez           elementwise rhs vectors
%          xyz           vertex coordinate vector
%          ev            element mapping matrix
%          ebound3D      element face boundary matrix
%          xl_m,yl_m,zl_m,xl_s,yl_s,zl_s  element coordinates 
%   output
%          elerr_p       element error estimate with BC correction
%
%  calls functions gausspoints_oned, gausspoints_threed deriv3D
%  IFISS function: GP 9 June 2022
% Copyright (c) 2022  G.Papanikos, C.E. Powell, D.J. Silvester
ngpt=2;
[oneg,onew] = gausspoints_oned(ngpt);
[s,t,l,wt] = gausspoints_threed(oneg,onew);
nngpt=ngpt^3;

x=xyz(:,1); y=xyz(:,2); z= xyz(:,3);
nel=length(ev(:,1));
lev=[ev,ev(:,end)]; elerr_p=elerror;
aez= zeros(nel,19,19);
ael =zeros(nel,8,8,8);
% inner-inner loop over subdivided elements
for subelt=1:8
    
    xl_m = permute(xl_s(:,subelt,:),[1,3,2]);
    yl_m = permute(yl_s(:,subelt,:),[1,3,2]);
    zl_m = permute(zl_s(:,subelt,:),[1,3,2]);
    % loop over Gauss points
    for igpt = 1:nngpt
        sigpt=s(igpt);
        tigpt=t(igpt);
        ligpt=l(igpt);
        wght=wt(igpt);
        % evaluate derivatives etc
        [~,invjac_v,~,dphidx_v,dphidy_v,dphidz_v] = deriv3D(sigpt,tigpt,ligpt,xl_m,yl_m,zl_m);
        ael(:,:,:,subelt) =ael(:,:,:,subelt)+ wght*(dphidx_v.*permute(dphidx_v,[1 3 2]) + .....
            dphidy_v.*permute(dphidy_v,[1 3 2]) + .....
            dphidz_v.*permute(dphidz_v,[1 3 2])).*invjac_v(:);
        
        % end of Gauss point loop
    end
 
end
ael  = permute(ael,[1,4,2,3]);
% manual assembly of subelement contributions
% --------------------------------------------------- first face ------------------------
% first edge,  node number 1                          % second edge,  node number 7             % third edge node number 15  
   aez(:,1,1)  = ael(:,1,2,2) + ael(:,2,1,1);       aez(:,7,7)   = ael(:,2,3,3) + ael(:,3,2,2);     aez(:,15,15) = ael(:,3,4,4) + ael(:,4,3,3);
   aez(:,1,2)  = ael(:,2,1,6);                      aez(:,7,1)   = ael(:,2,3,1);                    aez(:,15,13) = ael(:,4,3,1);
   aez(:,1,4)  = ael(:,1,2,5);                      aez(:,7,15)  = ael(:,3,2,4);                    aez(:,15,7)  = ael(:,3,4,2); 
   aez(:,1,7)  = ael(:,2,1,3);                      aez(:,7,8)   = ael(:,2,3,7) + ael(:,3,2,6);     aez(:,15,6)  = ael(:,3,4,1) + ael(:,4,3,2);
   aez(:,1,13) = ael(:,1,2,4);                      aez(:,7,5)   = ael(:,2,3,5);                    aez(:,15,16) = ael(:,3,4,7);
   aez(:,1,5)  = ael(:,1,2,6) + ael(:,2,1,5);       aez(:,7,16)  = ael(:,3,2,7);                    aez(:,15,18) = ael(:,4,3,8);
   aez(:,1,6)  = ael(:,1,2,3) + ael(:,2,1,4);       aez(:,7,6)   = ael(:,2,3,4) + ael(:,3,2,1);     aez(:,15,19) = ael(:,3,4,8) + ael(:,4,3,7);
   aez(:,1,14) = ael(:,1,2,7) + ael(:,2,1,8);       aez(:,7,14)  = ael(:,2,3,8) + ael(:,3,2,5);     aez(:,15,14) = ael(:,3,4,5) + ael(:,4,3,6);
   aez(:,1,12) = ael(:,1,2,8);                      aez(:,7,2)   = ael(:,2,3,6);                    aez(:,15,12) = ael(:,4,3,5);
   aez(:,1,8)  = ael(:,2,1,7);                      aez(:,7,19)  = ael(:,3,2,8);                    aez(:,15,8)  = ael(:,3,4,6);
                                     
% fourth edge, node 13                               %cedroid, node number 6 
  aez(:,13,13) = ael(:,4,1,1) + ael(:,1,4,4);       aez(:,6,6)   = ael(:,1,3,3) + ael(:,2,4,4) + ael(:,3,1,1) + ael(:,4,2,2);
  aez(:,13,1)  = ael(:,1,4,2);                      aez(:,6,1)   = ael(:,1,3,2) + ael(:,2,4,1);
  aez(:,13,15) = ael(:,4,1,3);                      aez(:,6,7)   = ael(:,2,4,3) + ael(:,3,1,2);
  aez(:,13,6)  = ael(:,4,1,2) + ael(:,1,4,3);       aez(:,6,15)  = ael(:,3,1,4) + ael(:,4,2,3);
  aez(:,13,4)  = ael(:,1,4,5);                      aez(:,6,13)  = ael(:,4,2,1) + ael(:,1,3,4);
  aez(:,13,18) = ael(:,4,1,8);                      aez(:,6,5)   = ael(:,1,3,6) + ael(:,2,4,5);
  aez(:,13,12) = ael(:,1,4,8) + ael(:,4,1,5);       aez(:,6,8)   = ael(:,2,4,7) + ael(:,3,1,6);
  aez(:,13,14) = ael(:,1,4,7) + ael(:,4,1,6);       aez(:,6,19)  = ael(:,3,1,8) + ael(:,4,2,7);
  aez(:,13,5)  = ael(:,1,4,6);                      aez(:,6,12)  = ael(:,4,2,5) + ael(:,1,3,8);
  aez(:,13,19) = ael(:,4,1,7);                      aez(:,6,14)  = ael(:,1,3,7) + ael(:,2,4,8) + ael(:,3,1,5) + ael(:,4,2,6);
                                                    aez(:,6,18)  = ael(:,4,2,8);
                                                    aez(:,6,4)   = ael(:,1,3,5);
                                                    aez(:,6,2)   = ael(:,2,4,6);
                                                    aez(:,6,16)  = ael(:,3,1,7);
                                                   
% -------------------------------------------------- Second face ------------------------                     
% first node number 7 is a common node of the faces 1 and 2 and is already
% assembled
   
% second node, number 16                             % third node, number 9                            % fourth node, number 2
  aez(:,16,16) = ael(:,3,7,7) + ael(:,7,3,3);       aez(:,9,9)  = ael(:,6,7,7) + ael(:,7,6,6);       aez(:,2,2)  = ael(:,2,6,6) + ael(:,6,2,2);         
  aez(:,16,9)  = ael(:,7,3,6);                      aez(:,9,2)  = ael(:,6,7,2);                      aez(:,2,9)  = ael(:,6,2,7);
  aez(:,16,7)  = ael(:,3,7,2);                      aez(:,9,16) = ael(:,7,6,3);                      aez(:,2,7)  = ael(:,2,6,3);
  aez(:,16,8)  = ael(:,3,7,6)+ ael(:,7,3,2);        aez(:,9,8)  = ael(:,6,7,3) + ael(:,7,6,2);       aez(:,2,8)  = ael(:,2,6,7) + ael(:,6,2,3);
  aez(:,16,15) = ael(:,3,7,4);                      aez(:,9,3)  = ael(:,6,7,5);                      aez(:,2,1)  = ael(:,2,6,1);
  aez(:,16,17) = ael(:,7,3,8);                      aez(:,9,17) = ael(:,7,6,8);                      aez(:,2,3)  = ael(:,6,2,5);
  aez(:,16,19) = ael(:,3,7,8) + ael(:,7,3,4);       aez(:,9,5)  = ael(:,6,7,1);                      aez(:,2,5)  = ael(:,2,6,5) + ael(:,6,2,1);
  aez(:,16,14) = ael(:,3,7,5) + ael(:,7,3,1);       aez(:,9,19) = ael(:,7,6,4);                      aez(:,2,14) = ael(:,2,6,8) + ael(:,6,2,4);
  aez(:,16,10) = ael(:,7,3,5);                      aez(:,9,14) = ael(:,6,7,4) + ael(:,7,6,1);       aez(:,2,10) = ael(:,6,2,8);
  aez(:,16,6)  = ael(:,3,7,1);                      aez(:,9,10) = ael(:,6,7,8) + ael(:,7,6,5);       aez(:,2,6)  = ael(:,2,6,4); 
 
% cedroid, node number 8   
  aez(:,8,8)  = ael(:,2,7,7) + ael(:,6,3,3) + ael(:,3,6,6) + ael(:,7,2,2);
  aez(:,8,2)  = ael(:,2,7,6) + ael(:,6,3,2);
  aez(:,8,16) = ael(:,3,6,7) + ael(:,7,2,3);
  aez(:,8,9)  = ael(:,6,3,7) + ael(:,7,2,6);
  aez(:,8,7)  = ael(:,2,7,3) + ael(:,3,6,2);
  aez(:,8,14) = ael(:,2,7,8) + ael(:,3,6,5) + ael(:,6,3,4) + ael(:,7,2,1);
  aez(:,8,19) = ael(:,3,6,8) + ael(:,7,2,4);
  aez(:,8,5)  = ael(:,2,7,5) + ael(:,6,3,1);
  aez(:,8,3)  = ael(:,6,3,5);
  aez(:,8,10) = ael(:,6,3,8) + ael(:,7,2,5);
  aez(:,8,17) = ael(:,7,2,8);
  aez(:,8,1)  = ael(:,2,7,1);
  aez(:,8,6)  = ael(:,2,7,4)+ ael(:,3,6,1);
  aez(:,8,15) = ael(:,3,6,4);
 
% -------------------------------------------------- third face ------------------------
% second node number 9 is a common node between face 2 and 3 and is already
% assembled

% first node, number 3                               % third node 17,                             %  fourth node number 11
  aez(:,3,3)   = ael(:,5,6,6) + ael(:,6,5,5);        aez(:,17,17) = ael(:,7,8,8) + ael(:,8,7,7);     aez(:,11,11) = ael(:,5,8,8) + ael(:,8,5,5);        
  aez(:,3,11)  = ael(:,5,6,8);                       aez(:,17,11) = ael(:,8,7,5);                    aez(:,11,17) = ael(:,8,5,7);
  aez(:,3,9)   = ael(:,6,5,7);                       aez(:,17,9)  = ael(:,7,8,6);                    aez(:,11,3)  = ael(:,5,8,6);
  aez(:,3,10)  = ael(:,5,6,7) + ael(:,6,5,8);        aez(:,17,10) = ael(:,7,8,5) + ael(:,8,7,6);     aez(:,11,10) = ael(:,5,8,7) + ael(:,8,5,6);
  aez(:,3,4)   = ael(:,5,6,1);                       aez(:,17,18) = ael(:,8,7,4);                    aez(:,11,4)  = ael(:,5,8,1);
  aez(:,3,2)   = ael(:,6,5,2);                       aez(:,17,19) = ael(:,7,8,4)+ael(:,8,7,3);       aez(:,11,12) = ael(:,5,8,4) + ael(:,8,5,1);
  aez(:,3,5)   = ael(:,5,6,2)+ael(:,6,5,1);          aez(:,17,16) = ael(:,7,8,3);                    aez(:,11,18) = ael(:,8,5,4); 
  aez(:,3,12)  = ael(:,5,6,4);                       aez(:,17,12) = ael(:,8,7,1);                    aez(:,11,5)  = ael(:,5,8,2);
  aez(:,3,14)  = ael(:,5,6,3)+ael(:,6,5,4);          aez(:,17,14) = ael(:,7,8,1) + ael(:,8,7,2);     aez(:,11,14) = ael(:,5,8,3)+ ael(:,8,5,2);
  aez(:,3,8)   = ael(:,6,5,3);                       aez(:,17,8)  = ael(:,7,8,2);                    aez(:,11,19) = ael(:,8,5,3);

% centroid node number 10
  aez(:,10,10) = ael(:,5,7,7) + ael(:,6,8,8) + ael(:,7,5,5) + ael(:,8,6,6);
  aez(:,10,3)  = ael(:,5,7,6) + ael(:,6,8,5);
  aez(:,10,9)  = ael(:,6,8,7) + ael(:,7,5,6);
  aez(:,10,17) = ael(:,7,5,8) + ael(:,8,6,7);
  aez(:,10,11) = ael(:,8,6,5) + ael(:,5,7,8);
  aez(:,10,5)  = ael(:,5,7,2) + ael(:,6,8,1);                                % found missing element connection a6_81
  aez(:,10,14) = ael(:,5,7,3) + ael(:,6,8,4) + ael(:,7,5,1) + ael(:,8,6,2);
  aez(:,10,19) = ael(:,8,6,3) + ael(:,7,5,4);                                % found missing element connection a7,54
  aez(:,10,4)  = ael(:,5,7,1);
  aez(:,10,12) = ael(:,5,7,4) + ael(:,8,6,1);
  aez(:,10,18) = ael(:,8,6,4);
  aez(:,10,2)  = ael(:,6,8,2);
  aez(:,10,8)  = ael(:,6,8,3)+ael(:,7,5,2); 
  aez(:,10,16) = ael(:,7,5,3);
  
  % -------------------------------------------------- fourth face ------------------------
  % node 13 is common with face 1 and node 11 with face 3 and are already
  % assembled
  
  %second node, number 18                            %fourth node, number 4                     % centroid node, number 12
  aez(:,18,18) = ael(:,4,8,8) + ael(:,8,4,4);       aez(:,4,4)  = ael(:,1,5,5) + ael(:,5,1,1);      aez(:,12,12) = ael(:,1,8,8) + ael(:,4,5,5) + ael(:,5,4,4) + ael(:,8,1,1);
  aez(:,18,11) = ael(:,8,4,5);                      aez(:,4,11) = ael(:,5,1,8);                     aez(:,12,13) = ael(:,1,8,4) + ael(:,4,5,1);
  aez(:,18,12) = ael(:,4,8,5) + ael(:,8,4,1);       aez(:,4,12) = ael(:,1,5,8) + ael(:,5,1,4);      aez(:,12,18) = ael(:,4,5,8) + ael(:,8,1,4);
  aez(:,18,13) = ael(:,4,8,1);                      aez(:,4,13) = ael(:,1,5,4);                     aez(:,12,11) = ael(:,5,4,8) + ael(:,8,1,5);
  aez(:,18,15) = ael(:,4,8,3);                      aez(:,4,3)  = ael(:,5,1,6);                     aez(:,12,4)  = ael(:,1,8,5) + ael(:,5,4,1);
  aez(:,18,19) = ael(:,4,8,7) + ael(:,8,4,3);       aez(:,4,5)  = ael(:,1,5,6) + ael(:,5,1,2);      aez(:,12,1)  = ael(:,1,8,2);
  aez(:,18,17) = ael(:,8,4,7);                      aez(:,4,1)  = ael(:,1,5,2);                     aez(:,12,6)  = ael(:,1,8,3) + ael(:,4,5,2);
  aez(:,18,6)  = ael(:,4,8,2);                      aez(:,4,6)  = ael(:,1,5,3);                     aez(:,12,15) = ael(:,4,5,3);
  aez(:,18,14) = ael(:,4,8,6) + ael(:,8,4,2);       aez(:,4,14) = ael(:,1,5,7) + ael(:,5,1,3);      aez(:,12,5)  = ael(:,1,8,6) + ael(:,5,4,2);  % add missing connection a5_42
  aez(:,18,10) = ael(:,8,4,6);                      aez(:,4,10) = ael(:,5,1,7);                     aez(:,12,14) = ael(:,1,8,7) + ael(:,4,5,6)+ael(:,5,4,3)+ael(:,8,1,2);
                                                                                                    aez(:,12,19) = ael(:,4,5,7) + ael(:,8,1,3);
                                                                                                    aez(:,12,3)  = ael(:,5,4,6);
                                                                                                    aez(:,12,10) = ael(:,5,4,7) + ael(:,8,1,6);
                                                                                                    aez(:,12,17) = ael(:,8,1,7);
                                                                                              
 % center nodes of fith sixth face and the centre node all 8 elements
 
 % node number 5 (center node of the fifth face),   node 19 (center node
 % of the sixth face) node 14 center node
 % ----------------------------------------------------------------------
 aez(:,5,5)  = ael(:,1,6,6) + ael(:,2,5,5) + ael(:,5,2,2) + ael(:,6,1,1);   aez(:,19,19) = ael(:,3,8,8) + ael(:,4,7,7) + ael(:,7,4,4) + ael(:,8,3,3);   aez(:,14,14) = ael(:,1,7,7) + ael(:,2,8,8) + ael(:,3,5,5) + ael(:,4,6,6) + ael(:,5,3,3) + ael(:,6,4,4) + ael(:,7,1,1)+ael(:,8,2,2); 
 aez(:,5,1)  = ael(:,1,6,2) + ael(:,2,5,1);                                 aez(:,19,15) = ael(:,3,8,4) + ael(:,4,7,3);                                 aez(:,14,6)  = ael(:,1,7,3) + ael(:,2,8,4) + ael(:,3,5,1) + ael(:,4,6,2);
 aez(:,5,2)  = ael(:,2,5,6) + ael(:,6,1,2);                                 aez(:,19,16) = ael(:,3,8,7) + ael(:,7,4,3);                                 aez(:,14,7)  = ael(:,2,8,3) + ael(:,3,5,2);
 aez(:,5,3)  = ael(:,5,2,6) + ael(:,6,1,5);                                 aez(:,19,17) = ael(:,7,4,8) + ael(:,8,3,7);                                 aez(:,14,8)  = ael(:,2,8,7) + ael(:,3,5,6) + ael(:,6,4,3) + ael(:,7,1,2);
 aez(:,5,4)  = ael(:,1,6,5) + ael(:,5,2,1);                                 aez(:,19,18) = ael(:,4,7,8) + ael(:,8,3,4);                                 aez(:,14,9)  = ael(:,6,4,7) + ael(:,7,1,6);
 aez(:,5,13) = ael(:,1,6,4);                                                aez(:,19,13) = ael(:,4,7,1);                                                aez(:,14,10) = ael(:,5,3,7) + ael(:,6,4,8) + ael(:,7,1,5) + ael(:,8,2,6);
 aez(:,5,6)  = ael(:,1,6,3) + ael(:,2,5,4);                                 aez(:,19,6)  = ael(:,4,7,2) + ael(:,3,8,1);                                 aez(:,14,11) = ael(:,5,3,8) + ael(:,8,2,5);
 aez(:,5,7)  = ael(:,2,5,3);                                                aez(:,19,7)  = ael(:,3,8,2);                                                aez(:,14,12) = ael(:,1,7,8) + ael(:,4,6,5) + ael(:,5,3,4) + ael(:,8,2,1);
 aez(:,5,12) = ael(:,1,6,8) + ael(:,5,2,4);                                 aez(:,19,12) = ael(:,4,7,5) + ael(:,8,3,1);                                 aez(:,14,13) = ael(:,1,7,4) + ael(:,4,6,1);
 aez(:,5,14) = ael(:,1,6,7) + ael(:,2,5,8) + ael(:,5,2,3) + ael(:,6,1,4);   aez(:,19,14) = ael(:,3,8,5) + ael(:,4,7,6) + ael(:,7,4,1) + ael(:,8,3,2);   aez(:,14,15) = ael(:,3,5,4) + ael(:,4,6,3);
 aez(:,5,8)  = ael(:,2,5,7) + ael(:,6,1,3);                                 aez(:,19,8)  = ael(:,3,8,6) + ael(:,7,4,2);                                 aez(:,14,19) = ael(:,3,5,8) + ael(:,4,6,7) + ael(:,7,1,4) + ael(8,2,3);
 aez(:,5,11) = ael(:,5,2,8);                                                aez(:,19,11) = ael(:,8,3,5);                                                aez(:,14,17) = ael(:,7,1,8) + ael(:,8,2,7);
 aez(:,5,10) = ael(:,5,2,7)+ael(:,6,1,8);                                   aez(:,19,10) = ael(:,8,3,6)+ael(:,7,4,5);                                   aez(:,14,3)  = ael(:,5,3,6) + ael(:,6,4,5);
 aez(:,5,9)  = ael(:,6,1,7);                                                aez(:,19,9)  = ael(:,7,4,6);                                                aez(:,14,5)  = ael(:,1,7,6) + ael(:,2,8,5) + ael(:,5,3,2) + ael(:,6,4,1);
                                                                                                                                                        aez(:,14,1)  = ael(:,1,7,2) + ael(:,2,8,1);                                                                                                                                                   
                                                                                                                                                        aez(:,14,2)  = ael(:,2,8,6) + ael(:,6,4,2);
                                                                                                                                                        aez(:,14,4)  = ael(:,1,7,5) + ael(:,5,3,1);  % added missing connection 14 - 4
                                                                                                                                                        aez(:,14,16) = ael(:,3,5,7) + ael(:,7,1,3);  % added missing connection 14 - 16
                                                                                                                                                        aez(:,14,18) = ael(:,4,6,8) + ael(:,8,2,4);  % added missing connection 14 - 18
                                                                                                                                             


nn=19;                                                                                                                                                       
mm=19;                                                                                                                                                        

tic
% recompute contributions from elements with Dirichlet boundaries
nbdf=length(ebound3D(:,1));
fbdy = zeros(nel,1);
face = zeros(nel,1);
% isolate boundary elements
for el = 1:nbdf
    ee = ebound3D(el,1);
    fbdy(ee) = fbdy(ee)+1; face(ee)=ebound3D(el,2);
end

% three face elements
k3=find(fbdy==3);
nel3b=length(k3);
% loop over two edge elements
for el = 1:nel3b
    el3e=k3(el);
    kk=find(ebound3D(:,1) == el3e);
    faces=ebound3D(kk,2);
    % set up original matrix and RHS vector
    ae=squeeze(aez(el3e,1:nn,1:mm));
    fe=fez(el3e,:)';
    % set up local coordinates and impose interpolated error as Dirichlet bc
    xl=x(lev(el3e,:)); yl=y(lev(el3e,:)); zl=z(lev(el3e,:));
    [bae,fe] = localbcfull_p3D(ae,fe,faces,xl,yl,zl);
    % solve local problem
    err=bae\fe;
    elerr_p(el3e,1) = err'*fe;
 
end
% end of element loop
%
%
% two face elements
k2=find(fbdy==2);
nel2b=length(k2);
% loop over two face elements
for el = 1:nel2b
    el2e=k2(el);
    kk=find(ebound3D(:,1) == el2e);
    faces=ebound3D(kk,2);
    % set up original matrix and RHS vector
    ae=squeeze(aez(el2e,1:nn,1:mm));
    fe=fez(el2e,:)';
    % set up local coordinates and impose interpolated error as Dirichlet bc
    xl=x(lev(el2e,:)); yl=y(lev(el2e,:)); zl=z(lev(el2e,:));
    [bae,fe] = localbcfull_p3D(ae,fe,faces,xl,yl,zl);
    % solve local problem
    err=bae\fe;
    elerr_p(el2e,1) = err'*fe;
end
% end of element loop
%
% one edge elements
k1=find(fbdy==1);
nel1b=length(k1);
% loop over one face elements
for el = 1:nel1b
    el1e=k1(el);
    kk=find(ebound3D(:,1) == el1e);
    faces=ebound3D(kk,2);
    % set up original matrix and RHS vector
    fe=fez(el1e,:)';
    ae=squeeze(aez(el1e,1:nn,1:mm));
    % set up local coordinates and impose interpolated error as Dirichlet bc
    xl=x(lev(el1e,:)); yl=y(lev(el1e,:)); zl=z(lev(el1e,:));
    [bae,fe] = localbcfull_p3D(ae,fe,faces,xl,yl,zl);
    % solve local problem
    err=bae\fe;
    elerr_p(el1e,1) = err'*fe;
end
% end of element loop
%
etime=toc;
fprintf('error boundary correction took %6.3e seconds\n',etime)
return