function macrogridplot3D(xyz,mv,bound3D,mbound3D)
%MACROGRIDPLOT3D hexahedral macroelement grid verification
%   macrogridplot3D(xyz,mv,bound3D,mbound3D);
% IFISS function: GP; 9 June 2022.
% Copyright (c)  2022  G.Papanikos,  C.E. Powell, D.J. Silvester
fprintf('\nSubdivision logistics ..\n')
nvtx=length(xyz(:,1));
fprintf('  %g nodes \n',nvtx)
nelement=length(mv(:,1));
fprintf('  %g (nxnxn macro)elements \n',nelement)
nboundvtx=length(bound3D);
fprintf('  %g nodes on Dirichlet boundary \n',nboundvtx)
nboundedge=length(mbound3D(:,1));
fprintf('  %g macroelement faces on Dirichlet boundary\n\n',nboundedge)

adj=sparse(nvtx,nvtx); adx=sparse(nvtx,nvtx);
for i=1:length(mv(:,1))
    
    adj(mv(i,1),mv(i,2))=1;
    adj(mv(i,1),mv(i,9))=1;
    adj(mv(i,1),mv(i,12))=1;
    adj(mv(i,1),mv(i,21))=1;
    adj(mv(i,1),mv(i,4))=1;
    adj(mv(i,1),mv(i,5))=1;
    
    
    adj(mv(i,2),mv(i,3))=1;     adx(mv(i,20),mv(i,19))=1;
    adj(mv(i,2),mv(i,1))=1;      adx(mv(i,20),mv(i,21))=1;
    adj(mv(i,2),mv(i,9))=1;      adx(mv(i,20),mv(i,12))=1;
    adj(mv(i,2),mv(i,6))=1;      adx(mv(i,20),mv(i,26))=1;
    adj(mv(i,2),mv(i,10))=1;    adx(mv(i,20),mv(i,22))=1;
    adj(mv(i,2),mv(i,15))=1;     adx(mv(i,20),mv(i,16))=1;
    
    adj(mv(i,3),mv(i,4))=1;     adx(mv(i,14),mv(i,15))=1;
    adj(mv(i,3),mv(i,15))=1;    adx(mv(i,14),mv(i,9))=1;
    adj(mv(i,3),mv(i,23))=1;    adx(mv(i,14),mv(i,23))=1;
    adj(mv(i,3),mv(i,2))=1;     adx(mv(i,14),mv(i,21))=1;
    adj(mv(i,3),mv(i,24))=1;    adx(mv(i,14),mv(i,22))=1;
    adj(mv(i,3),mv(i,7))=1;     adx(mv(i,14),mv(i,18))=1;
    
    adj(mv(i,4),mv(i,1))=1;     adx(mv(i,16),mv(i,24))=1; adx(mv(i,18),mv(i,25))=1;
    adj(mv(i,4),mv(i,3))=1;     adx(mv(i,16),mv(i,10))=1; adx(mv(i,18),mv(i,11))=1;
    adj(mv(i,4),mv(i,21))=1;    adx(mv(i,16),mv(i,17))=1; adx(mv(i,18),mv(i,19))=1;
    adj(mv(i,4),mv(i,23))=1;    adx(mv(i,16),mv(i,15))=1; adx(mv(i,18),mv(i,17))=1;
    adj(mv(i,4),mv(i,26))=1;    adx(mv(i,16),mv(i,22))=1; adx(mv(i,18),mv(i,22))=1;
    adj(mv(i,4),mv(i,8))=1;     adx(mv(i,16),mv(i,20))=1; adx(mv(i,18),mv(i,14))=1;
    
    adj(mv(i,5),mv(i,6))=1;    adx(mv(i,13),mv(i,22))=1;  adx(mv(i,22),mv(i,13))=1; adx(mv(i,27),mv(i,22))=1;
    adj(mv(i,5),mv(i,11))=1;   adx(mv(i,13),mv(i,27))=1;  adx(mv(i,22),mv(i,16))=1; adx(mv(i,27),mv(i,26))=1;
    adj(mv(i,5),mv(i,1))=1;    adx(mv(i,13),mv(i,10))=1;  adx(mv(i,22),mv(i,18))=1; adx(mv(i,27),mv(i,25))=1;
    adj(mv(i,5),mv(i,12))=1;   adx(mv(i,13),mv(i,12))=1;  adx(mv(i,22),mv(i,27))=1; adx(mv(i,27),mv(i,23))=1;
    adj(mv(i,5),mv(i,19))=1;   adx(mv(i,13),mv(i,9))=1;   adx(mv(i,22),mv(i,20))=1; adx(mv(i,27),mv(i,24))=1;
    adj(mv(i,5),mv(i,8))=1;    adx(mv(i,13),mv(i,11))=1;  adx(mv(i,22),mv(i,14))=1; adx(mv(i,27),mv(i,13))=1;
    
    
    adj(mv(i,6),mv(i,7))=1;
    adj(mv(i,6),mv(i,11))=1;
    adj(mv(i,6),mv(i,2))=1;
    adj(mv(i,6),mv(i,5))=1;
    adj(mv(i,6),mv(i,17))=1;
    adj(mv(i,6),mv(i,10))=1;
    
    adj(mv(i,7),mv(i,8))=1;
    adj(mv(i,7),mv(i,3))=1;
    adj(mv(i,7),mv(i,25))=1;
    adj(mv(i,7),mv(i,24))=1;
    adj(mv(i,7),mv(i,17))=1;
    adj(mv(i,7),mv(i,6))=1;
    
    
    adj(mv(i,8),mv(i,5))=1;
    adj(mv(i,8),mv(i,4))=1;
    adj(mv(i,8),mv(i,25))=1;
    adj(mv(i,8),mv(i,26))=1;
    adj(mv(i,8),mv(i,19))=1;
    adj(mv(i,8),mv(i,7))=1;
end
%% define element faces
adjb=sparse(nvtx,nvtx);
%bottom boundary face z=-1


if isempty(find(mbound3D(:,2)==7,1))
    %% Cube domain
    k1=find(mbound3D(:,2)==1)';
    for k=mbound3D(k1)
        
        adjb(mv(k,1),mv(k,9))=1;
        adjb(mv(k,9),mv(k,2))=1;
        adjb(mv(k,2),mv(k,15))=1;
        adjb(mv(k,15),mv(k,3))=1;
        adjb(mv(k,3),mv(k,23))=1;
        adjb(mv(k,23),mv(k,4))=1;
        adjb(mv(k,4),mv(k,21))=1;
        adjb(mv(k,21),mv(k,1))=1;
        
        adjb(mv(k,15),mv(k,14))=1;
        adjb(mv(k,14),mv(k,21))=1;
        adjb(mv(k,14),mv(k,9))=1;
        adjb(mv(k,14),mv(k,23))=1;
    end
    
    % x=1
    % right boundary face
    k2=find(mbound3D(:,2)==2)';
    for k=mbound3D(k2)
        
        adjb(mv(k,2),mv(k,15))=1;
        adjb(mv(k,15),mv(k,3))=1;
        adjb(mv(k,3),mv(k,24))=1;
        adjb(mv(k,24),mv(k,7))=1;
        adjb(mv(k,7),mv(k,17))=1;
        adjb(mv(k,17),mv(k,16))=1;
        adjb(mv(k,16),mv(k,15))=1;
        adjb(mv(k,17),mv(k,6))=1;
        adjb(mv(k,6),mv(k,10))=1;
        adjb(mv(k,6),mv(k,2))=1;
        
    end
    %
    % top boundary face z=1
    k3=find(mbound3D(:,2)==3)';
    for k=mbound3D(k3)
        adjb(mv(k,5),mv(k,11))=1;
        adjb(mv(k,11),mv(k,6))=1;
        adjb(mv(k,6),mv(k,17))=1;
        adjb(mv(k,17),mv(k,7))=1;
        adjb(mv(k,7),mv(k,25))=1;
        adjb(mv(k,25),mv(k,8))=1;
        adjb(mv(k,8),mv(k,19))=1;
        adjb(mv(k,19),mv(k,5))=1;
        
        adjb(mv(k,25),mv(k,18))=1;
        adjb(mv(k,18),mv(k,11))=1;
        adjb(mv(k,18),mv(k,17))=1;
        adjb(mv(k,19),mv(k,18))=1;
    end
    % x=-1
    % back boundary face
    
    k4=find(mbound3D(:,2)==4)';
    for k=mbound3D(k4)
        adjb(mv(k,1),mv(k,4))=1;
        adjb(mv(k,4),mv(k,8))=1;
        adjb(mv(k,8),mv(k,5))=1;
        adjb(mv(k,5),mv(k,1))=1;
        
        adjb(mv(k,8),mv(k,19))=1;
        adjb(mv(k,19),mv(k,20))=1;
        adjb(mv(k,20),mv(k,12))=1;
        adjb(mv(k,20),mv(k,26))=1;
        adjb(mv(k,20),mv(k,21))=1;
    end
    
    
    % left boundary face y=-1
    k5=find(mbound3D(:,2)==5)';
    for k=mbound3D(k5)
        
        adjb(mv(k,1),mv(k,9))=1;
        adjb(mv(k,9),mv(k,2))=1;
        adjb(mv(k,9),mv(k,13))=1;
        
        adjb(mv(k,2),mv(k,10))=1;
        adjb(mv(k,10),mv(k,6))=1;
        adjb(mv(k,10),mv(k,13))=1;
        
        adjb(mv(k,6),mv(k,11))=1;
        adjb(mv(k,11),mv(k,13))=1;
        adjb(mv(k,11),mv(k,5))=1;
        adjb(mv(k,5),mv(k,12))=1;
        adjb(mv(k,12),mv(k,13))=1;
        adjb(mv(k,12),mv(k,1))=1;
        
    end
    
    %y=1
    % back boundary face
    k6=find(mbound3D(:,2)==6)';
    for k=mbound3D(k6)
        adjb(mv(k,4),mv(k,23))=1;
        adjb(mv(k,23),mv(k,3))=1;
        adjb(mv(k,3),mv(k,24))=1;
        adjb(mv(k,24),mv(k,7))=1;
        adjb(mv(k,7),mv(k,25))=1;
        adjb(mv(k,25),mv(k,8))=1;
        adjb(mv(k,8),mv(k,26))=1;
        adjb(mv(k,26),mv(k,4))=1;
        
        adjb(mv(k,26),mv(k,27))=1;
        adjb(mv(k,27),mv(k,24))=1;
        adjb(mv(k,27),mv(k,25))=1;
        adjb(mv(k,27),mv(k,23))=1;
        
    end
    
else
    %% stair domain
    % z=0
    k1=find(mbound3D(:,2)==1)';
    for k=mbound3D(k1)
        
        adjb(mv(k,1),mv(k,9))=1;
        adjb(mv(k,9),mv(k,2))=1;
        adjb(mv(k,2),mv(k,15))=1;
        adjb(mv(k,15),mv(k,3))=1;
        adjb(mv(k,3),mv(k,23))=1;
        adjb(mv(k,23),mv(k,4))=1;
        adjb(mv(k,4),mv(k,21))=1;
        adjb(mv(k,21),mv(k,1))=1;
        
        adjb(mv(k,15),mv(k,14))=1;
        adjb(mv(k,14),mv(k,21))=1;
        adjb(mv(k,14),mv(k,9))=1;
        adjb(mv(k,14),mv(k,23))=1;
    end
    % x=0
    k2=find(mbound3D(:,2)==2)';
    for k=mbound3D(k2)
        adjb(mv(k,1),mv(k,4))=1;
        adjb(mv(k,4),mv(k,8))=1;
        adjb(mv(k,8),mv(k,5))=1;
        adjb(mv(k,5),mv(k,1))=1;
        
        adjb(mv(k,8),mv(k,19))=1;
        adjb(mv(k,19),mv(k,20))=1;
        adjb(mv(k,20),mv(k,12))=1;
        adjb(mv(k,20),mv(k,26))=1;
        adjb(mv(k,20),mv(k,21))=1;
    end
    %  z=-1
    k3=find(mbound3D(:,2)==3)';
    for k=mbound3D(k3)
        
        adjb(mv(k,1),mv(k,9))=1;
        adjb(mv(k,9),mv(k,2))=1;
        adjb(mv(k,2),mv(k,15))=1;
        adjb(mv(k,15),mv(k,3))=1;
        adjb(mv(k,3),mv(k,23))=1;
        adjb(mv(k,23),mv(k,4))=1;
        adjb(mv(k,4),mv(k,21))=1;
        adjb(mv(k,21),mv(k,1))=1;
        
        adjb(mv(k,15),mv(k,14))=1;
        adjb(mv(k,14),mv(k,21))=1;
        adjb(mv(k,14),mv(k,9))=1;
        adjb(mv(k,14),mv(k,23))=1;
    end
    
    % x=1
    % right boundary face
    k4=find(mbound3D(:,2)==4)';
    for k=mbound3D(k4)
        
        adjb(mv(k,2),mv(k,15))=1;
        adjb(mv(k,15),mv(k,3))=1;
        adjb(mv(k,3),mv(k,24))=1;
        adjb(mv(k,24),mv(k,7))=1;
        adjb(mv(k,7),mv(k,17))=1;
        adjb(mv(k,17),mv(k,16))=1;
        adjb(mv(k,16),mv(k,15))=1;
        adjb(mv(k,17),mv(k,6))=1;
        adjb(mv(k,6),mv(k,10))=1;
        adjb(mv(k,6),mv(k,2))=1;
        
    end
    
    % top boundary face z=1
    k5=find(mbound3D(:,2)==5)';
    for k=mbound3D(k5)
        adjb(mv(k,5),mv(k,11))=1;
        adjb(mv(k,11),mv(k,6))=1;
        adjb(mv(k,6),mv(k,17))=1;
        adjb(mv(k,17),mv(k,7))=1;
        adjb(mv(k,7),mv(k,25))=1;
        adjb(mv(k,25),mv(k,8))=1;
        adjb(mv(k,8),mv(k,19))=1;
        adjb(mv(k,19),mv(k,5))=1;
        
        adjb(mv(k,25),mv(k,18))=1;
        adjb(mv(k,18),mv(k,11))=1;
        adjb(mv(k,18),mv(k,17))=1;
        adjb(mv(k,19),mv(k,18))=1;
    end
    
    
    % x=-1
    % back boundary face
    
    k6=find(mbound3D(:,2)==6)';
    for k=mbound3D(k6)
        adjb(mv(k,1),mv(k,4))=1;
        adjb(mv(k,4),mv(k,8))=1;
        adjb(mv(k,8),mv(k,5))=1;
        adjb(mv(k,5),mv(k,1))=1;
        
        adjb(mv(k,8),mv(k,19))=1;
        adjb(mv(k,19),mv(k,20))=1;
        adjb(mv(k,20),mv(k,12))=1;
        adjb(mv(k,20),mv(k,26))=1;
        adjb(mv(k,20),mv(k,21))=1;
    end
    
    
    % left boundary face y=-1
    k7=find(mbound3D(:,2)==7)';
    for k=mbound3D(k7)
        
        adjb(mv(k,1),mv(k,9))=1;
        adjb(mv(k,9),mv(k,2))=1;
        adjb(mv(k,9),mv(k,13))=1;
        
        adjb(mv(k,2),mv(k,10))=1;
        adjb(mv(k,10),mv(k,6))=1;
        adjb(mv(k,10),mv(k,13))=1;
        
        adjb(mv(k,6),mv(k,11))=1;
        adjb(mv(k,11),mv(k,13))=1;
        adjb(mv(k,11),mv(k,5))=1;
        adjb(mv(k,5),mv(k,12))=1;
        adjb(mv(k,12),mv(k,13))=1;
        adjb(mv(k,12),mv(k,1))=1;
        
    end
    
    %y=1
    % back boundary face
    k8=find(mbound3D(:,2)==8)';
    for k=mbound3D(k8)
        adjb(mv(k,4),mv(k,23))=1;
        adjb(mv(k,23),mv(k,3))=1;
        adjb(mv(k,3),mv(k,24))=1;
        adjb(mv(k,24),mv(k,7))=1;
        adjb(mv(k,7),mv(k,25))=1;
        adjb(mv(k,25),mv(k,8))=1;
        adjb(mv(k,8),mv(k,26))=1;
        adjb(mv(k,26),mv(k,4))=1;
        
        adjb(mv(k,26),mv(k,27))=1;
        adjb(mv(k,27),mv(k,24))=1;
        adjb(mv(k,27),mv(k,25))=1;
        adjb(mv(k,27),mv(k,23))=1;
        
    end
    
    
    
    
end

figure(1)
gplot3(adj,xyz,'b');
hold on
gplot3(adx,xyz,'c');
%stnode=int2str([1:nvtx]');
%text(xyz(:,1),xyz(:,2),xyz(:,3),stnode)
axis('equal'),axis('off')

xlabel('x')
ylabel('y')
zlabel('z')

title('Indices of nodes of the macroelement grid')
hold off



figure(2)
gplot3(adj,xyz,'b');
hold on
gplot3(adx,xyz,'c');
gplot3(adjb,xyz,'r');

xlabel('x')
ylabel('y')
zlabel('z')

xyzbd=xyz(bound3D,:);
stbd=int2str([1:nboundvtx]');
text(xyzbd(:,1),xyzbd(:,2),xyzbd(:,3),stbd,'color','black')
title('Indices of nodes on the Dirichlet boundary')
axis('equal'),axis('off')
hold off



end
