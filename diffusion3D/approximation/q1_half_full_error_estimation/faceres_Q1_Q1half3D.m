function [faceres] = faceres_Q1_Q1half3D(xyz,ev,ebound3D,q1sol3D,fcx,hx,hy,hz)
%FACERES_Q1_Q1half3D  computes Q1(h/2) face residuals for Q1 solution
%  [faceres] = faceres_Q1_Q1half3D(xyz,ev,ebound3D,q1sol3D,fcx,hx,hy,hz)
%   input
%          xyz           vertex coordinate vector
%          ev            element mapping matrix
%          ebound3D      element face boundary matrix
%          q1sol3D       stochastic Q1 solution vector
%          fcx           element face connectivity array
%          hx,hy,hz      element mesh sizes
%   output
%          faceres      edge residuals
%
%   calls functions q1fluxjmps3D_modified, gausspoints_oned, stoch_gauss_coeff, deriv3D qderiv
%   SIFISS function: GP; 09 June 2022
% Copyright (c) 2022 G. Papanikos, C.E. Powell, D.J. Silvester

x=xyz(:,1); y=xyz(:,2); z=xyz(:,3); nvtx=length(x);
nel=length(ev(:,1));

% construct the integration rule
ngpt=2;
[oneg,onew] = gausspoints_oned(ngpt);
[s,t,wt] = gausspoints_twod(oneg,onew);
nngpt=ngpt^2;
% initialisation
faceres = zeros(nel,19);
xl_v=nan(nel,8); yl_v=nan(nel,8); zl_v=nan(nel,8);
 
% inner loop over elements
xl_v = x(ev);
yl_v = y(ev);
zl_v = z(ev);

% loop over Gauss points
for igpt=1:nngpt
    sigpt=s(igpt);
    sigpt_l=(sigpt-1.0e0)/2; sigpt_r=(sigpt+1.0e0)/2;
    tigpt=t(igpt);
    tigpt_l=(tigpt-1.0e0)/2; tigpt_r=(tigpt+1.0e0)/2;
    
    wigpt=wt(igpt);
    
    % face 1
    [~,~,phi_v_1,~,~,~] = deriv3D(sigpt,tigpt,-1,xl_v,yl_v,zl_v);
    % face 2
    
    [~,~,phi_v_2,~,~,~] = deriv3D(1,sigpt,tigpt,xl_v,yl_v,zl_v);
  
    % face 3
    [~,~,phi_v_3,~,~,~] = deriv3D(sigpt,tigpt,1,xl_v,yl_v,zl_v);
    
    %face 4
    
    [~,~,phi_v_4,~,~,~] = deriv3D(-1,sigpt,tigpt,xl_v,yl_v,zl_v);
    
    %face 5
    
    [~,~,phi_v_5,~,~,~]  = deriv3D(sigpt,-1,tigpt,xl_v,yl_v,zl_v);
    %face 6
    
    [~,~,phi_v_6,~,~,~]  = deriv3D(sigpt,1,tigpt,xl_v,yl_v,zl_v);
    
    
    jmp_ll = q1fluxjmps3D_modified(q1sol3D,fcx,xyz,ev,ebound3D,sigpt_l,tigpt_l);
    jmp_rr = q1fluxjmps3D_modified(q1sol3D,fcx,xyz,ev,ebound3D,sigpt_r,tigpt_r);
    jmp_lr = q1fluxjmps3D_modified(q1sol3D,fcx,xyz,ev,ebound3D,sigpt_l,tigpt_r);
    jmp_rl = q1fluxjmps3D_modified(q1sol3D,fcx,xyz,ev,ebound3D,sigpt_r,tigpt_l);
    %----------------------------------------------------------
    %------------------node numbering -------------------------
    % 9 10 11 12 13  14  15 16 17 18 19 20 21 22 23 24 25 26 27
    % 1  2  3  4  5   6   7  8  9 10 11 12 13 14 15 16 17 18 19
    
    
    % -------------- face 1 -----------------------------------------------
    faceres(:,1)  = faceres(:,1)  +  (1/2)*(1/2)*wigpt*jmp_ll(:,1).*phi_v_1(:,2).*hx(:).*hy(:)./8;
    faceres(:,1)  = faceres(:,1)  +  (1/2)*(1/2)*wigpt*jmp_rl(:,1).*phi_v_1(:,1).*hx(:).*hy(:)./8;
    
    faceres(:,7)  = faceres(:,7)  +  (1/2)*(1/2)*wigpt*jmp_rl(:,1).*phi_v_1(:,3).*hx(:).*hy(:)./8;
    faceres(:,7)  = faceres(:,7)  +  (1/2)*(1/2)*wigpt*jmp_rr(:,1).*phi_v_1(:,2).*hx(:).*hy(:)./8;
    
    faceres(:,15) = faceres(:,15) +  (1/2)*(1/2)*wigpt*jmp_rr(:,1).*phi_v_1(:,4).*hx(:).*hy(:)./8;
    faceres(:,15) = faceres(:,15) +  (1/2)*(1/2)*wigpt*jmp_lr(:,1).*phi_v_1(:,3).*hx(:).*hy(:)./8;
    
    faceres(:,13) = faceres(:,13) +  (1/2)*(1/2)*wigpt*jmp_lr(:,1).*phi_v_1(:,1).*hx(:).*hy(:)./8;
    faceres(:,13) = faceres(:,13) +  (1/2)*(1/2)*wigpt*jmp_ll(:,1).*phi_v_1(:,4).*hx(:).*hy(:)./8;
    
    faceres(:,6)  = faceres(:,6)  +  (1/2)*(1/2)*wigpt*jmp_ll(:,1).*phi_v_1(:,3).*hx(:).*hy(:)./8;
    faceres(:,6)  = faceres(:,6)  +  (1/2)*(1/2)*wigpt*jmp_rl(:,1).*phi_v_1(:,4).*hx(:).*hy(:)./8;
    faceres(:,6)  = faceres(:,6)  +  (1/2)*(1/2)*wigpt*jmp_rr(:,1).*phi_v_1(:,1).*hx(:).*hy(:)./8;
    faceres(:,6)  = faceres(:,6)  +  (1/2)*(1/2)*wigpt*jmp_lr(:,1).*phi_v_1(:,2).*hx(:).*hy(:)./8;
    
    % ---------------- face 2 ---------------------------------------------
   
    
    faceres(:,7)  = faceres(:,7)  +  (1/2)*(1/2)*wigpt*jmp_ll(:,2).*phi_v_2(:,3).*hy(:).*hz(:)./8;
    faceres(:,7)  = faceres(:,7)  +  (1/2)*(1/2)*wigpt*jmp_rl(:,2).*phi_v_2(:,2).*hy(:).*hz(:)./8;
    
    faceres(:,16) = faceres(:,16) +  (1/2)*(1/2)*wigpt*jmp_rl(:,2).*phi_v_2(:,7).*hy(:).*hz(:)./8;
    faceres(:,16) = faceres(:,16) +  (1/2)*(1/2)*wigpt*jmp_rr(:,2).*phi_v_2(:,3).*hy(:).*hz(:)./8;
    
     
    faceres(:,9)  = faceres(:,9)  +  (1/2)*(1/2)*wigpt*jmp_rr(:,2).*phi_v_2(:,6).*hy(:).*hz(:)./8;
    faceres(:,9)  = faceres(:,9)  +  (1/2)*(1/2)*wigpt*jmp_lr(:,2).*phi_v_2(:,7).*hy(:).*hz(:)./8;
    
    
    faceres(:,2)  = faceres(:,2)  +  (1/2)*(1/2)*wigpt*jmp_lr(:,2).*phi_v_2(:,2).*hy(:).*hz(:)./8;
    faceres(:,2)  = faceres(:,2)  +  (1/2)*(1/2)*wigpt*jmp_ll(:,2).*phi_v_2(:,6).*hy(:).*hz(:)./8;
    
  
    faceres(:,8)  = faceres(:,8)  +  (1/2)*(1/2)*wigpt*jmp_ll(:,2).*phi_v_2(:,7).*hy(:).*hz(:)./8;
    faceres(:,8)  = faceres(:,8)  +  (1/2)*(1/2)*wigpt*jmp_rl(:,2).*phi_v_2(:,6).*hy(:).*hz(:)./8;
    faceres(:,8)  = faceres(:,8)  +  (1/2)*(1/2)*wigpt*jmp_rr(:,2).*phi_v_2(:,2).*hy(:).*hz(:)./8;
    faceres(:,8)  = faceres(:,8)  +  (1/2)*(1/2)*wigpt*jmp_lr(:,2).*phi_v_2(:,3).*hy(:).*hz(:)./8;
    
    % --------------- face 3 ----------------------------------------------
    faceres(:,3)  = faceres(:,3)  +  (1/2)*(1/2)*wigpt*jmp_ll(:,3).*phi_v_3(:,6).*hx(:).*hy(:)./8;
    faceres(:,3)  = faceres(:,3)  +  (1/2)*(1/2)*wigpt*jmp_rl(:,3).*phi_v_3(:,5).*hx(:).*hy(:)./8;
    
    faceres(:,9)  = faceres(:,9)  +  (1/2)*(1/2)*wigpt*jmp_rl(:,3).*phi_v_3(:,7).*hx(:).*hy(:)./8;
    faceres(:,9)  = faceres(:,9)  +  (1/2)*(1/2)*wigpt*jmp_rr(:,3).*phi_v_3(:,6).*hx(:).*hy(:)./8;
    
    faceres(:,17) = faceres(:,17) +  (1/2)*(1/2)*wigpt*jmp_rr(:,3).*phi_v_3(:,8).*hx(:).*hy(:)./8;
    faceres(:,17) = faceres(:,17) +  (1/2)*(1/2)*wigpt*jmp_lr(:,3).*phi_v_3(:,7).*hx(:).*hy(:)./8;
    
    faceres(:,11) = faceres(:,11) +  (1/2)*(1/2)*wigpt*jmp_lr(:,3).*phi_v_3(:,5).*hx(:).*hy(:)./8;
    faceres(:,11) = faceres(:,11) +  (1/2)*(1/2)*wigpt*jmp_ll(:,3).*phi_v_3(:,8).*hx(:).*hy(:)./8;
    
    faceres(:,10) = faceres(:,10) +  (1/2)*(1/2)*wigpt*jmp_ll(:,3).*phi_v_3(:,7).*hx(:).*hy(:)./8;
    faceres(:,10) = faceres(:,10) +  (1/2)*(1/2)*wigpt*jmp_rl(:,3).*phi_v_3(:,8).*hx(:).*hy(:)./8;
    faceres(:,10) = faceres(:,10) +  (1/2)*(1/2)*wigpt*jmp_rr(:,3).*phi_v_3(:,5).*hx(:).*hy(:)./8;
    faceres(:,10) = faceres(:,10) +  (1/2)*(1/2)*wigpt*jmp_lr(:,3).*phi_v_3(:,6).*hx(:).*hy(:)./8;
    
    %  ---------- face 4 --------------------------------------------------
    faceres(:,13)  = faceres(:,13) +  (1/2)*(1/2)*wigpt*jmp_ll(:,4).*phi_v_4(:,4).*hy(:).*hz(:)./8;
    faceres(:,13)  = faceres(:,13) +  (1/2)*(1/2)*wigpt*jmp_rl(:,4).*phi_v_4(:,1).*hy(:).*hz(:)./8;
    
    faceres(:,18)  = faceres(:,18) +  (1/2)*(1/2)*wigpt*jmp_rl(:,4).*phi_v_4(:,8).*hy(:).*hz(:)./8;
    faceres(:,18)  = faceres(:,18) +  (1/2)*(1/2)*wigpt*jmp_rr(:,4).*phi_v_4(:,4).*hy(:).*hz(:)./8;
    
    faceres(:,11)  = faceres(:,11) +  (1/2)*(1/2)*wigpt*jmp_rr(:,4).*phi_v_4(:,5).*hy(:).*hz(:)./8;
    faceres(:,11)  = faceres(:,11) +  (1/2)*(1/2)*wigpt*jmp_lr(:,4).*phi_v_4(:,8).*hy(:).*hz(:)./8;
    
    faceres(:,4)   = faceres(:,4)  +  (1/2)*(1/2)*wigpt*jmp_lr(:,4).*phi_v_4(:,1).*hy(:).*hz(:)./8;
    faceres(:,4)   = faceres(:,4)  +  (1/2)*(1/2)*wigpt*jmp_ll(:,4).*phi_v_4(:,5).*hy(:).*hz(:)./8;
    
    faceres(:,12)  = faceres(:,12) +  (1/2)*(1/2)*wigpt*jmp_ll(:,4).*phi_v_4(:,8).*hy(:).*hz(:)./8;
    faceres(:,12)  = faceres(:,12) +  (1/2)*(1/2)*wigpt*jmp_rl(:,4).*phi_v_4(:,5).*hy(:).*hz(:)./8;
    faceres(:,12)  = faceres(:,12) +  (1/2)*(1/2)*wigpt*jmp_rr(:,4).*phi_v_4(:,1).*hy(:).*hz(:)./8;
    faceres(:,12)  = faceres(:,12) +  (1/2)*(1/2)*wigpt*jmp_lr(:,4).*phi_v_4(:,4).*hy(:).*hz(:)./8;
    
    % --------- face 5 ----------------------------------------------------
    faceres(:,1)  = faceres(:,1) +   (1/2)*(1/2)*wigpt*jmp_ll(:,5).*phi_v_5(:,2).*hx(:).*hz(:)./8;
    faceres(:,1)  = faceres(:,1) +   (1/2)*(1/2)*wigpt*jmp_rl(:,5).*phi_v_5(:,1).*hx(:).*hz(:)./8;
    
    faceres(:,2)  = faceres(:,2) +   (1/2)*(1/2)*wigpt*jmp_rl(:,5).*phi_v_5(:,6).*hx(:).*hz(:)./8;
    faceres(:,2)  = faceres(:,2) +   (1/2)*(1/2)*wigpt*jmp_rr(:,5).*phi_v_5(:,2).*hx(:).*hz(:)./8;
    
    faceres(:,3)  = faceres(:,3) +   (1/2)*(1/2)*wigpt*jmp_rr(:,5).*phi_v_5(:,5).*hx(:).*hz(:)./8;
    faceres(:,3)  = faceres(:,3) +   (1/2)*(1/2)*wigpt*jmp_lr(:,5).*phi_v_5(:,6).*hx(:).*hz(:)./8;
    
    faceres(:,4)  = faceres(:,4) +   (1/2)*(1/2)*wigpt*jmp_lr(:,5).*phi_v_5(:,1).*hx(:).*hz(:)./8;
    faceres(:,4)  = faceres(:,4) +   (1/2)*(1/2)*wigpt*jmp_ll(:,5).*phi_v_5(:,5).*hx(:).*hz(:)./8;
    
    faceres(:,5)  = faceres(:,5) +   (1/2)*(1/2)*wigpt*jmp_ll(:,5).*phi_v_5(:,6).*hx(:).*hz(:)./8;
    faceres(:,5)  = faceres(:,5) +   (1/2)*(1/2)*wigpt*jmp_rl(:,5).*phi_v_5(:,5).*hx(:).*hz(:)./8;
    faceres(:,5)  = faceres(:,5) +   (1/2)*(1/2)*wigpt*jmp_rr(:,5).*phi_v_5(:,1).*hx(:).*hz(:)./8;
    faceres(:,5)  = faceres(:,5) +   (1/2)*(1/2)*wigpt*jmp_lr(:,5).*phi_v_5(:,2).*hx(:).*hz(:)./8;
    
    % ------- face 6 ------------------------------------------------------
    faceres(:,15)  = faceres(:,15) +   (1/2)*(1/2)*wigpt*jmp_ll(:,6).*phi_v_6(:,3).*hx(:).*hz(:)./8;
    faceres(:,15)  = faceres(:,15) +   (1/2)*(1/2)*wigpt*jmp_rl(:,6).*phi_v_6(:,4).*hx(:).*hz(:)./8;
    
    faceres(:,16)  = faceres(:,16) +   (1/2)*(1/2)*wigpt*jmp_rl(:,6).*phi_v_6(:,7).*hx(:).*hz(:)./8;
    faceres(:,16)  = faceres(:,16) +   (1/2)*(1/2)*wigpt*jmp_rr(:,6).*phi_v_6(:,3).*hx(:).*hz(:)./8;
    
    faceres(:,17)  = faceres(:,17) +   (1/2)*(1/2)*wigpt*jmp_rr(:,6).*phi_v_6(:,8).*hx(:).*hz(:)./8;
    faceres(:,17)  = faceres(:,17) +   (1/2)*(1/2)*wigpt*jmp_lr(:,6).*phi_v_6(:,7).*hx(:).*hz(:)./8;
    
    faceres(:,18)  = faceres(:,18) +   (1/2)*(1/2)*wigpt*jmp_lr(:,6).*phi_v_6(:,4).*hx(:).*hz(:)./8;
    faceres(:,18)  = faceres(:,18) +   (1/2)*(1/2)*wigpt*jmp_ll(:,6).*phi_v_6(:,8).*hx(:).*hz(:)./8;
    
    faceres(:,19)  = faceres(:,19) +   (1/2)*(1/2)*wigpt*jmp_ll(:,6).*phi_v_6(:,7).*hx(:).*hz(:)./8;
    faceres(:,19)  = faceres(:,19) +   (1/2)*(1/2)*wigpt*jmp_rl(:,6).*phi_v_6(:,8).*hx(:).*hz(:)./8;
    faceres(:,19)  = faceres(:,19) +   (1/2)*(1/2)*wigpt*jmp_rr(:,6).*phi_v_6(:,4).*hx(:).*hz(:)./8;
    faceres(:,19)  = faceres(:,19) +   (1/2)*(1/2)*wigpt*jmp_lr(:,6).*phi_v_6(:,3).*hx(:).*hz(:)./8;
end           % end of Gauss points loop

return