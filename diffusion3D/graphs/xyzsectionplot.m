function xyzsectionplot(domain,qmethod,xyz,x_it,yyref,zzref,fig,mkr)
%XYZSECTIONPLOT plots scalar solution along line in three dimensions
%   xyzsectionplot(domain,qmethod,xyz,x_it,yyref,zzref,fig,'.b');
%   input
%          domain     domain identifier 
%          qmethod    approximation method 
%          xyz        nodal coordinate vector
%          x_it       solution vector
%          yyref      y-location of grid line
%          zzref      z-location of grid line
%          fig        figure number
%          mkr        plotting character
%
%   IFISS function: DJS; 9 August 2022: 22 September 2022
% Copyright (c) 2022 G.Papanikos, C.E. Powell, D.J. Silvester

nvtx=length(xyz(:,1));
if domain==3 %-------- xyz left-hand coordinate for vertical bore hole
   fprintf('x-section analysis | y = %8.4e  z = %8.4e\n',yyref,zzref)
   kkall=find(abs((xyz(:,2))-yyref)+abs((xyz(:,3))-zzref)<10*eps);
else %------- xzy left-hand coordinate
fprintf('x-section analysis | z = %8.4e  y = %8.4e\n',zzref,yyref)
   kkall=find(abs((xyz(:,2))-zzref)+abs((xyz(:,3))-yyref)<10*eps);
end
lkk=length(kkall);
kk=kkall; %kk=kkall(lkk-npts+1:lkk);
xref=xyz(kk,1); yref=xyz(kk,2); zref=xyz(kk,3);
if isempty(kk), error('Oops.. check location of grid line'), end
uxref=x_it(kk);
figure(fig)
plot(xref,uxref,mkr), axis('square'), title('x-section solution')
[xref,uxref]
return
