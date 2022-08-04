function [bae,fe] = localbcfull_p3D(ae,fe,faces,xl,yl,zl)
%LOCALBCFULL_P3D imposes Dirichlet BC for Poisson error estimator for full set
%of bubble functions   
%  [bae,fe] = localbcfull_p3D(ae,fe,edges,xl,yl,zl);
%  inputs:
%          ae           Poisson problem matrix
%          fe           rhs vector
%          faces        boundary face vector 
%          xl, yl, zl   vertex coordinates  
%  outputs:
%          bae       Poisson problem matrix
%          fe        rhs vector
%
%  calls function: specific_bc3D
%  IFISS function: GP; 09 June 2022
% Copyright (c)  2022 G.Papanikos, C.E. Powell, D.J. Silvester
bae=ae; ffe=fe;
nvtx = length(fe); nbd=length(faces);
zero_col=zeros(nvtx,5); zero_row=zeros(5,nvtx);
error = zeros(size(fe));
boundary_nodes = [1 7 15 13 6; 7 16 9 2 8; 3 9 17 11 10; 13 18 11 4 12; 1 2 3 4 5; 15 16 17 18 19];
%
%% loop over boundary edges
for bd=1:nbd
    ek=faces(bd);
% compute boundary edge coordinates
    if ek==1
       zbd(1:9)= min(zl);  
       
       xbd(1)= min(xl);              ybd(1)= min(yl);  
       xbd(3)= max(xl);              ybd(3)= min(yl);
       xbd(2)=0.5*(xbd(1)+xbd(3));   ybd(2)=0.5*(ybd(1)+ybd(3));  
       xbd(5)= max(xl);              ybd(5)= max(yl);      
       xbd(4)=0.5*(xbd(5)+xbd(3));   ybd(4)=0.5*(ybd(5)+ybd(3)); 
       xbd(7)= min(xl);              ybd(7)= max(yl);             
       xbd(6)=0.5*(xbd(7)+xbd(5));   ybd(6)=0.5*(ybd(7)+ybd(5));
       xbd(8)=0.5*(xbd(7)+xbd(1));   ybd(8)=0.5*(ybd(7)+ybd(1)); 
       xbd(9)=0.5*(xbd(2)+xbd(6));   ybd(9)=0.5*(ybd(2)+ybd(6)); 
    
    elseif ek ==2     
       xbd(1:9)= max(xl);  
       
       ybd(1)= min(yl);              zbd(1)= min(zl);  
       ybd(3)= max(yl);              zbd(3)= min(zl);
       ybd(2)=0.5*(ybd(1)+ybd(3));   zbd(2)=0.5*(zbd(1)+zbd(3));  
       ybd(5)= max(yl);              zbd(5)= max(zl);      
       ybd(4)=0.5*(ybd(5)+ybd(3));   zbd(4)=0.5*(zbd(5)+zbd(3)); 
       ybd(7)= min(yl);              zbd(7)= max(zl);             
       ybd(6)=0.5*(ybd(7)+ybd(5));   zbd(6)=0.5*(zbd(7)+zbd(5));
       ybd(8)=0.5*(ybd(7)+ybd(1));   zbd(8)=0.5*(zbd(7)+zbd(1));
       ybd(9)=0.5*(ybd(2)+ybd(6));   zbd(9)=0.5*(zbd(2)+zbd(6));  
    elseif ek ==3
       zbd(1:9)= max(zl);  
       
       xbd(1)= min(xl);              ybd(1)= min(yl);  
       xbd(3)= max(xl);              ybd(3)= min(yl);
       xbd(2)=0.5*(xbd(1)+xbd(3));   ybd(2)=0.5*(ybd(1)+ybd(3));  
       xbd(5)= max(xl);              ybd(5)= max(yl);      
       xbd(4)=0.5*(xbd(5)+xbd(3));   ybd(4)=0.5*(ybd(5)+ybd(3)); 
       xbd(7)= min(xl);              ybd(7)= max(yl);             
       xbd(6)=0.5*(xbd(7)+xbd(5));   ybd(6)=0.5*(ybd(7)+ybd(5));
       xbd(8)=0.5*(xbd(7)+xbd(1));   ybd(8)=0.5*(ybd(7)+ybd(1)); 
       xbd(9)=0.5*(xbd(2)+xbd(6));   ybd(9)=0.5*(ybd(2)+ybd(6));
       
    elseif ek ==4
       xbd(1:9)= min(xl);  
       
       ybd(1)= min(yl);              zbd(1)= min(zl);  
       ybd(3)= max(yl);              zbd(3)= min(zl);
       ybd(2)=0.5*(ybd(1)+ybd(3));   zbd(2)=0.5*(zbd(1)+zbd(3));  
       ybd(5)= max(yl);              zbd(5)= max(zl);      
       ybd(4)=0.5*(ybd(5)+ybd(3));   zbd(4)=0.5*(zbd(5)+zbd(3)); 
       ybd(7)= min(yl);              zbd(7)= max(zl);             
       ybd(6)=0.5*(ybd(7)+ybd(5));   zbd(6)=0.5*(zbd(7)+zbd(5));
       ybd(8)=0.5*(ybd(7)+ybd(1));   zbd(8)=0.5*(zbd(7)+zbd(1));   
       ybd(9)=0.5*(ybd(2)+ybd(6));   zbd(9)=0.5*(zbd(2)+zbd(6));
       
       
    elseif ek ==5
       ybd(1:9)= min(yl);  
       
       xbd(1)= min(xl);              zbd(1)= min(zl);  
       xbd(3)= max(xl);              zbd(3)= min(zl);
       xbd(2)=0.5*(xbd(1)+xbd(3));   zbd(2)=0.5*(zbd(1)+zbd(3));  
       xbd(5)= max(xl);              zbd(5)= max(zl);      
       xbd(4)=0.5*(xbd(5)+xbd(3));   zbd(4)=0.5*(zbd(5)+zbd(3)); 
       xbd(7)= min(xl);              zbd(7)= max(zl);             
       xbd(6)=0.5*(xbd(7)+xbd(5));   zbd(6)=0.5*(zbd(7)+zbd(5));
       xbd(8)=0.5*(xbd(7)+xbd(1));   zbd(8)=0.5*(zbd(7)+zbd(1)); 
       xbd(9)=0.5*(xbd(2)+xbd(6));   zbd(9)=0.5*(zbd(2)+zbd(6));
    else
       ybd(1:9)= max(yl);  
       
       xbd(1)= min(xl);              zbd(1)= min(zl);  
       xbd(3)= max(xl);              zbd(3)= min(zl);
       xbd(2)=0.5*(xbd(1)+xbd(3));   zbd(2)=0.5*(zbd(1)+zbd(3));  
       xbd(5)= max(xl);              zbd(5)= max(zl);      
       xbd(4)=0.5*(xbd(5)+xbd(3));   zbd(4)=0.5*(zbd(5)+zbd(3)); 
       xbd(7)= min(xl);              zbd(7)= max(zl);             
       xbd(6)=0.5*(xbd(7)+xbd(5));   zbd(6)=0.5*(zbd(7)+zbd(5));
       xbd(8)=0.5*(xbd(7)+xbd(1));   zbd(8)=0.5*(zbd(7)+zbd(1));
       xbd(9)=0.5*(xbd(2)+xbd(6));   zbd(9)=0.5*(zbd(2)+zbd(6));
    end
    
  % compute interpolated boundary error at 5 nodes  
    
    bd_nodes = boundary_nodes(ek,:);
  
    bc=specific_bc3D(xbd,ybd,zbd);
    error(bd_nodes(1))  =  bc(2)-0.5*(bc(1)+bc(3));
    error(bd_nodes(2))  =  bc(4)-0.5*(bc(3)+bc(5));
    error(bd_nodes(3))  =  bc(6)-0.5*(bc(5)+bc(7));
    error(bd_nodes(4))  =  bc(8)-0.5*(bc(1)+bc(7));
    error(bd_nodes(5))  =  bc(9)-1/4*(bc(1)+bc(3)+bc(5)+bc(7));
    
%% impose boundary condition without modifying the other equations
%% DJS/DK mod
    bae(:,bd_nodes)=zero_col; bae(bd_nodes,:)=zero_row;   
    bae(bd_nodes(1),bd_nodes(1))=1;  
    bae(bd_nodes(2),bd_nodes(2))=1;
    bae(bd_nodes(3),bd_nodes(3))=1;
    bae(bd_nodes(4),bd_nodes(4))=1;
    bae(bd_nodes(5),bd_nodes(5))=1;
    fe(bd_nodes)=error(bd_nodes) ;  
end
return
