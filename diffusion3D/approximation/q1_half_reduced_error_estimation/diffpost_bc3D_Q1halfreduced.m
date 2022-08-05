function [elerr_p] = diffpost_bc3D_Q1halfreduced(elerror,fez,xyz,ev,ebound3D,xl_m,yl_m,zl_m,xl_s,yl_s,zl_s)
%DIFFPOST_Q1HALFREDUCED postprocesses local Poisson error estimator 
%   [elerr_p] = diffpost_bc3D_Q1halfreduced(elerror,fez,xyz,ev,ebound3D,xl_m,yl_m,zl_m,xl_s,yl_s,zl_s) 
%   inputs:
%          elerror  	 element error estimate (without BC imposed)         
%          fez           elementwise rhs vectors
%          xyz           vertex coordinate vector  
%          ev            element mapping matrix
%          ebound3D      element face boundary matrix 
%	   xl_m,yl_m,zl_m,xl_s,yl_s,zl_s element coordinates
%   output
%          elerr_p      element error estimate with BC correction
%
%  Calls functions gausspoints_oned, gausspoints_threed
%  Rebuilt in IFISS3.1 to retain compatibility with earlier versions
%  IFISS function: GP 2021;
%  Copyright (c)   G.Papanikos, C.E. Powell, D.J. Silvester

     ngpt=2;
      [oneg,onew] = gausspoints_oned(ngpt);
      [s,t,l,wt] = gausspoints_threed(oneg,onew);
      nngpt=ngpt^3; 

      x=xyz(:,1); y=xyz(:,2); z= xyz(:,3);
      nel=length(ev(:,1));
      lev=[ev,ev(:,end)]; elerr_p=elerror;
      aez= zeros(nel,7,7);
      ael =zeros(nel,8,8,8);
% inner-inner loop over subdivided elements
  for subelt=1:8
      
          xl_m = permute(xl_s(:,subelt,:),[1 3 2]);
          yl_m = permute(yl_s(:,subelt,:),[1 3 2]);
          zl_m = permute(zl_s(:,subelt,:),[1 3 2]);

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
% end of subdivided element loop
  end
  ael  = permute(ael,[1,4,2,3]);
% manual assembly of subelement contributions
% ------------------------------- first face ------------------------

%cedroid, node number 1
aez(:,1,1)  = ael(:,1,3,3) + ael(:,2,4,4) + ael(:,3,1,1) + ael(:,4,2,2);
aez(:,1,2)  = ael(:,2,4,7) + ael(:,3,1,6);
aez(:,1,4)  = ael(:,4,2,5) + ael(:,1,3,8);
aez(:,1,5)  = ael(:,1,3,6) + ael(:,2,4,5);
aez(:,1,6)  = ael(:,3,1,8) + ael(:,4,2,7);
aez(:,1,7)  = ael(:,1,3,7) + ael(:,2,4,8) + ael(:,3,1,5) + ael(:,4,2,6);  %added missing connection

% -------------------------------- Second face ------------------------
% cedroid, node number 2
aez(:,2,1)  = ael(:,2,7,4) + ael(:,3,6,1);
aez(:,2,2)  = ael(:,2,7,7) + ael(:,6,3,3) + ael(:,3,6,6) + ael(:,7,2,2);
aez(:,2,3)  = ael(:,6,3,8) + ael(:,7,2,5);
aez(:,2,5)  = ael(:,2,7,5) + ael(:,6,3,1);
aez(:,2,6)  = ael(:,3,6,8) + ael(:,7,2,4);
aez(:,2,7)  = ael(:,2,7,8) + ael(:,3,6,5) + ael(:,6,3,4) + ael(:,7,2,1);

% -------------------------------- third face ------------------------

% centroid node number 3
aez(:,3,2) = ael(:,6,8,3) + ael(:,7,5,2);
aez(:,3,3) = ael(:,5,7,7) + ael(:,6,8,8) + ael(:,7,5,5) + ael(:,8,6,6);
aez(:,3,4) = ael(:,5,7,4) + ael(:,8,6,1);
aez(:,3,5) = ael(:,5,7,2) + ael(:,6,8,1);                                    % added missing connection
aez(:,3,6) = ael(:,8,6,3) + ael(:,7,5,4);                                    % added missing connection
aez(:,3,7) = ael(:,5,7,3) + ael(:,6,8,4) + ael(:,7,5,1) + ael(:,8,6,2);

% -------------------------------- fourth face ------------------------
% centroid node, number 4
aez(:,4,1) = ael(:,1,8,3) + ael(:,4,5,2);
aez(:,4,3) = ael(:,5,4,7) + ael(:,8,1,6);
aez(:,4,4) = ael(:,1,8,8) + ael(:,4,5,5) + ael(:,5,4,4) + ael(:,8,1,1);
aez(:,4,5) = ael(:,1,8,6) + ael(:,5,2,4);                                   %added missing connection
aez(:,4,6) = ael(:,4,5,7) + ael(:,8,1,3);
aez(:,4,7) = ael(:,1,8,7) + ael(:,4,5,6) + ael(:,5,4,3) + ael(:,8,1,2);

% -------------------------------fifth face --------------------------
% centroid node, number 5

aez(:,5,1) = ael(:,1,6,3) + ael(:,2,5,4);
aez(:,5,2) = ael(:,2,5,7) + ael(:,6,1,3);
aez(:,5,3) = ael(:,5,2,7) + ael(:,6,1,8);
aez(:,5,4) = ael(:,1,6,8) + ael(:,5,2,4);
aez(:,5,5) = ael(:,1,6,6) + ael(:,2,5,5) + ael(:,5,2,2) + ael(:,6,1,1);
aez(:,5,7) = ael(:,1,6,7) + ael(:,2,5,8) + ael(:,5,2,3) + ael(:,6,1,4);


% -------------------------------sixth face --------------------------
% centroid node, number 6

aez(:,6,1)  = ael(:,4,7,2) + ael(:,3,8,1);
aez(:,6,2)  = ael(:,3,8,6) + ael(:,7,4,2);
aez(:,6,3)  = ael(:,8,3,6) + ael(:,7,4,5);
aez(:,6,4)  = ael(:,4,7,5) + ael(:,8,3,1);
aez(:,6,6)  = ael(:,3,8,8) + ael(:,4,7,7) + ael(:,7,4,4) + ael(:,8,3,3);
aez(:,6,7)  = ael(:,3,8,5) + ael(:,4,7,6) + ael(:,7,4,1) + ael(:,8,3,2);

% ----------------  mid node  of all elements --------------------------- 
% number 7
aez(:,7,1)  = ael(:,1,7,3) + ael(:,2,8,4) + ael(:,3,5,1) + ael(:,4,6,2);
aez(:,7,2)  = ael(:,2,8,7) + ael(:,3,5,6) + ael(:,6,4,3) + ael(:,7,1,2);
aez(:,7,3)  = ael(:,5,3,7) + ael(:,6,4,8) + ael(:,7,1,5) + ael(:,8,2,6);
aez(:,7,4)  = ael(:,1,7,8) + ael(:,4,6,5) + ael(:,5,3,4) + ael(:,8,2,1);
aez(:,7,5)  = ael(:,1,7,6) + ael(:,2,8,5) + ael(:,5,3,2) + ael(:,6,4,1);
aez(:,7,6)  = ael(:,3,5,8) + ael(:,4,6,7) + ael(:,7,1,4) + ael(:,8,2,3);
aez(:,7,7)  = ael(:,1,7,7) + ael(:,2,8,8) + ael(:,3,5,5) + ael(:,4,6,6) + ael(:,5,3,3) + ael(:,6,4,4) + ael(:,7,1,1)+ael(:,8,2,2);                                                                                                                                                      aez(:,14,1)  = ael(:,1,7,2) + ael(:,2,8,1);                                                                                                                                                   
                                                                                                                                                        
                                                                                                                                                        
nn=7;                                                                                                                                                       
mm=7;                                                                                                                                                                                                                       
      
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
     [bae,fe] = localbc_p3D(ae,fe,faces,xl,yl,zl);
% solve local problem
      err=bae\fe;
      elerr_p(el3e,1) = err'*fe;
      
%       strct.faces{el} = [el3e,faces'];
%       strct.coord{el} =  [xl yl zl];
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
      [bae,fe] = localbc_p3D(ae,fe,faces,xl,yl,zl);
% solve local problem
      err=bae\fe;
      elerr_p(el2e,1) = err'*fe;
      end
% end of element loop
%
% one edge elements
      k1=find(fbdy==1);
      nel1b=length(k1);
% loop over one edge elements
      for el = 1:nel1b
      el1e=k1(el);
      kk=find(ebound3D(:,1) == el1e);
      faces=ebound3D(kk,2);
% set up original matrix and RHS vector 
      fe=fez(el1e,:)';
	  ae=squeeze(aez(el1e,1:nn,1:mm)); 
% set up local coordinates and impose interpolated error as Dirichlet bc
       xl=x(lev(el1e,:)); yl=y(lev(el1e,:)); zl=z(lev(el1e,:));
      [bae,fe] = localbc_p3D(ae,fe,faces,xl,yl,zl);
% solve local problem
      err=bae\fe;
      elerr_p(el1e,1) = err'*fe;
      end
% end of element loop
%
      etime=toc;
      fprintf('error boundary correction took %6.3e seconds\n',etime)  
 return
