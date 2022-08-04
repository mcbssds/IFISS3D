function [hx,hy,hz,efx,ebound3D] = facegen(xyz,ev,domain)
%FACEGEN face information for flux jump computation
%   [hx,hy,hz,efx] = facegen(xyz,ev,domain);
%   inputs:
%          xyz        nodal coordinate vector
%          ev         element mapping matrix
%          domain     cube, stair or borehole domain (1,2 or 3)
%   outputs:
%          hx,hy,hz   elementwise face areas
%          efx        element face index vector
%          ebound3D   boundary nodes
% IFISS function: GP 4th August 2022;
% Copyright (c)   G.Papanikos, C.E. Powell, D.J. Silvester

nvtx=length(xyz(:,1));
nel=length(ev(:,1));
x=xyz(:,1); y=xyz(:,2); z= xyz(:,3);

adj1 = sparse(nvtx,nvtx);
adj2 = sparse(nvtx,nvtx);
adj5 = sparse(nvtx,nvtx);


efc = zeros(nel,4,3);

%% compute edge to edge connection array ee
% initialise global matrices of faces
% evaluate element number on each face in turn
% and assemble into adjacency matrix


adj1  = adj1 + sparse(ev(:,1),ev(:,2),1:nel,nvtx,nvtx);
adj1  = adj1 + sparse(ev(:,2),ev(:,3),1:nel,nvtx,nvtx);
adj1  = adj1 + sparse(ev(:,3),ev(:,4),1:nel,nvtx,nvtx);
adj1  = adj1 + sparse(ev(:,4),ev(:,1),1:nel,nvtx,nvtx);

adj2  = adj2 + sparse(ev(:,2),ev(:,3),1:nel,nvtx,nvtx);
adj2  = adj2 + sparse(ev(:,3),ev(:,7),1:nel,nvtx,nvtx);
adj2  = adj2 + sparse(ev(:,7),ev(:,6),1:nel,nvtx,nvtx);
adj2  = adj2 + sparse(ev(:,6),ev(:,2),1:nel,nvtx,nvtx);

adj5  = adj5 + sparse(ev(:,1),ev(:,2),1:nel,nvtx,nvtx);
adj5  = adj5 + sparse(ev(:,2),ev(:,6),1:nel,nvtx,nvtx);
adj5  = adj5 + sparse(ev(:,6),ev(:,5),1:nel,nvtx,nvtx);
adj5  = adj5 + sparse(ev(:,5),ev(:,1),1:nel,nvtx,nvtx);


[pp1,kk1,ll1]=find(adj1);
A1 = [pp1 kk1 ll1];
[~,indx] = sort(A1(:,3),'ascend');
A1 = A1(indx,:);

[pp2,kk2,ll2]=find(adj2);
A2 = [pp2 kk2 ll2];
[~,indx] = sort(A2(:,3),'ascend');
A2 = A2(indx,:);

[pp5,kk5,ll5]=find(adj5);
A5 = [pp5 kk5 ll5];
[~,indx] = sort(A5(:,3),'ascend');
A5 = A5(indx,:);

elements = 1:nel;
ss  = [4*elements'-3, 4*elements'-2, 4*elements'-1, 4*elements'];
efc(:,:,1) = reshape(diag(adj1(A1(ss,2),A1(ss,1)))',nel,4);
efc(:,:,2) = reshape(diag(adj2(A2(ss,2),A2(ss,1)))',nel,4);
efc(:,:,3) = reshape(diag(adj5(A5(ss,2),A5(ss,1)))',nel,4);

efc = reshape(efc,[nel,12]);
if domain == 1 || domain == 3  % cube domain or borehole domain
    efc = efc(:,[1,10,8,4,5,6,7,3,9,2,11,12]);
    efc  = efc(:,[7,12,11,9,10,8]);
else  % stair domain
    efc = efc(:,[10,2,3,4,11,9,7,12,6,1,5,8]);
    efc = efc(:,7:12);
end

% inner loop over elements
xl_v = x(ev); yl_v = y(ev); zl_v = z(ev);
% compute local mesh sizes
hx=xl_v(:,2)-xl_v(:,1);
hy=yl_v(:,3)-yl_v(:,1);
hz=zl_v(:,5)-zl_v(:,1);

fprintf('checking face numbering and computing face areas ... ')
% centroid coordinates
for ielem=1:nel
    xc(ielem)=mean(x(ev(ielem,1:8))); yc(ielem)=mean(y(ev(ielem,1:8)));
    zc(ielem)=mean(z(ev(ielem,1:8)));
end
xyzp=[xc',yc',zc'];
%
% remove zero indices corresponding to boundary faces
ppk=[1:nel]'; onecol=ones(nel,1);
% find boundary faces
iie=efc==0;
efx=efc + ( efc==0).* [ppk,ppk,ppk,ppk,ppk,ppk];
%
k=1;
for i = 1:size(efx,1)
    k1  = find(efx(i,:)==i);
    sz = size(k1,2);
    e = i*ones(sz,1);
    if ~isempty(k1)
        ebound3Dnew{k,1} = [e, k1'];
        k=k+1;
    end
end
ebound3D= cell2mat(ebound3Dnew);


fprintf('done \n')
%
return
%--------------------------------------------------------------
function nn=mod6(n)
if n<7,
    nn=n;
else
    nn=mod(n,6);
end

%     end
% toc;
