function [err_sq_el,xx,fe,ae,xl_m,yl_m,zl_m,xl_s,yl_s,zl_s] = diffpost_Q1_Q1halfreduced3D(xyz,ev,ebound3D,q1sol3D,fcx,hx,hy,hz)
%DIFFPOST_Q1_Q1halfreduced3D computes Q1(h/2) reduced error estimator for Q1 solution
%employs elementwise reduced Q1-bubbles for midpoint centroid nodes
%   [zerr_sq_el,zelerr,zfe,zae] = stoch_diffpost_q1_ypx(xy,ev,ebound,x_gal,...
%                                  fcx,hx,hy,hz);
%   input
%          xyz          vertex coordinate vector
%          ev           element mapping matrix
%          ebound3D     element face boundary matrix
%          q1sol3D        Q1 solution vector
%          fcx          element face connectivity array
%          hx,hy,hz     element mesh sizes

%   output
%          err_sq_el   element error estimate
%          xx          elementwise error estimate
%          fe          elementwise rhs vectors
%          ae          LDLT factorized element matrices
%          xl_m,yl_m,zl_m,xl_s,yl_s,zl_s element coordinates
%
%   calls functions  gausspoints_oned, gausspoints_threed gauss_source3D
%   deriv3D
%   IFISS function: GP; 22 June 2022
% Copyright (c) 2022 G. Papanikos, C.E. Powell, D.J. Silvester

fprintf('computing Q1 error estimator...  \n')
x=xyz(:,1); y=xyz(:,2); z=xyz(:,3);
nel=length(ev(:,1));
xl_v=nan(nel,8); yl_v=nan(nel,8); zl_v=nan(nel,8);

%
% construct the integration rule
ngpt=3;
[oneg,onew] = gausspoints_oned(ngpt);
[s,t,l,wt] = gausspoints_threed(oneg,onew);
nngpt=ngpt^3;
t1=tic;
%
% inner loop over elements

xl_v = x(ev);
yl_v = y(ev);
zl_v = z(ev);


% initialise local stochastic and deterministic matrices
%  ae = zeros(nel,nn,nn);
ae = zeros(nel,7,7);
err_sq_el= zeros(nel,1);
elerr=zeros(7,nel);
ael = zeros(nel,8,8,8);
fde = zeros(nel,7);
fdem = zeros(nel,8,8);

xl_s = zeros(nel,8,8); yl_s = zeros(nel,8,8); zl_s = zeros(nel,8,8);
% compute local mid-edge coordinates

xface1(:,1)= 0.5*(xl_v(:,1)+xl_v(:,2));                                  yface1(:,1)= 0.5*(yl_v(:,1)+yl_v(:,2));                                   zface1(:,1) = zl_v(:,1);
xface1(:,2)= 0.5*(xl_v(:,2)+xl_v(:,3));                                  yface1(:,2)= 0.5*(yl_v(:,2)+yl_v(:,3));                                   zface1(:,2) = zl_v(:,1);
xface1(:,3)= 0.5*(xl_v(:,3)+xl_v(:,4));                                  yface1(:,3)= 0.5*(yl_v(:,3)+yl_v(:,4));                                   zface1(:,3) = zl_v(:,1);
xface1(:,4)= 0.5*(xl_v(:,1)+xl_v(:,4));                                  yface1(:,4)= 0.5*(yl_v(:,1)+yl_v(:,4));                                   zface1(:,4) = zl_v(:,1);
xcentrface1(:,1)= 0.25*( xface1(:,1)+ xface1(:,2)+ xface1(:,3)+ xface1(:,4)); ycentrface1(:,1)= 0.25*(yface1(:,1)+ yface1(:,2)+ yface1(:,3)+ yface1(:,4));   zcentrface1(:,1) = zl_v(:,1);

xface2(:,1)= xl_v(:,7);                                                  yface2(:,1)= 0.5*(yl_v(:,2)+yl_v(:,3));                                   zface2(:,1) = 0.5*(zl_v(:,2)+zl_v(:,3));
xface2(:,2)= xl_v(:,7);                                                  yface2(:,2)= 0.5*(yl_v(:,3)+yl_v(:,7));                                   zface2(:,2) = 0.5*(zl_v(:,3)+zl_v(:,7));
xface2(:,3)= xl_v(:,7);                                                  yface2(:,3)= 0.5*(yl_v(:,7)+yl_v(:,6));                                   zface2(:,3) = 0.5*(zl_v(:,7)+zl_v(:,6));
xface2(:,4)= xl_v(:,7);                                                  yface2(:,4)= 0.5*(yl_v(:,6)+yl_v(:,2));                                   zface2(:,4) = 0.5*(zl_v(:,6)+zl_v(:,2));
xcentrface2(:,1)= xl_v(:,7);                                             ycentrface2(:,1)= 0.25*(yface2(:,1)+ yface2(:,2)+ yface2(:,3)+ yface2(:,4));   zcentrface2(:,1) = 0.25*(zface2(:,1)+ zface2(:,2)+ zface2(:,3)+ zface2(:,4));

xface3(:,1)= 0.5*(xl_v(:,5)+xl_v(:,6));                                  yface3(:,1)= 0.5*(yl_v(:,5)+yl_v(:,6));                                   zface3(:,1) = zl_v(:,7);
xface3(:,2)= 0.5*(xl_v(:,6)+xl_v(:,7));                                  yface3(:,2)= 0.5*(yl_v(:,6)+yl_v(:,7));                                   zface3(:,2) = zl_v(:,7);
xface3(:,3)= 0.5*(xl_v(:,7)+xl_v(:,8));                                  yface3(:,3)= 0.5*(yl_v(:,7)+yl_v(:,8));                                   zface3(:,3) = zl_v(:,7);
xface3(:,4)= 0.5*(xl_v(:,8)+xl_v(:,5));                                  yface3(:,4)= 0.5*(yl_v(:,8)+yl_v(:,5));                                   zface3(:,4) = zl_v(:,7);
xcentrface3(:,1)= 0.25*( xface3(:,1)+ xface3(:,2)+ xface3(:,3)+ xface3(:,4)); ycentrface3(:,1)= 0.25*(yface3(:,1)+ yface3(:,2)+ yface3(:,3)+ yface3(:,4));   zcentrface3(:,1) = zl_v(:,7);

xface4(:,1)= xl_v(:,1);                                                  yface4(:,1)= 0.5*(yl_v(:,1)+yl_v(:,4));                                   zface4(:,1) = 0.5*(zl_v(:,1)+zl_v(:,4));
xface4(:,2)= xl_v(:,1);                                                  yface4(:,2)= 0.5*(yl_v(:,4)+yl_v(:,8));                                   zface4(:,2) = 0.5*(zl_v(:,4)+zl_v(:,8));
xface4(:,3)= xl_v(:,1);                                                  yface4(:,3)= 0.5*(yl_v(:,8)+yl_v(:,5));                                   zface4(:,3) = 0.5*(zl_v(:,8)+zl_v(:,5));
xface4(:,4)= xl_v(:,1);                                                  yface4(:,4)= 0.5*(yl_v(:,5)+yl_v(:,1));                                   zface4(:,4) = 0.5*(zl_v(:,5)+zl_v(:,1));
xcentrface4(:,1)= xl_v(:,1);                                             ycentrface4(:,1)= 0.25*(yface4(:,1)+ yface4(:,2)+ yface4(:,3)+ yface4(:,4));   zcentrface4(:,1) = 0.25*(zface4(:,1)+ zface4(:,2)+ zface4(:,3)+ zface4(:,4));

xface5(:,1)= 0.5*(xl_v(:,1)+xl_v(:,2));                                  yface5(:,1)= yl_v(:,1);                                                   zface5(:,1) = 0.5*(zl_v(:,1)+zl_v(:,2));
xface5(:,2)= 0.5*(xl_v(:,2)+xl_v(:,6));                                  yface5(:,2)= yl_v(:,1);                                                   zface5(:,2) = 0.5*(zl_v(:,2)+zl_v(:,6));
xface5(:,3)= 0.5*(xl_v(:,6)+xl_v(:,5));                                  yface5(:,3)= yl_v(:,1);                                                   zface5(:,3) = 0.5*(zl_v(:,6)+zl_v(:,5));
xface5(:,4)= 0.5*(xl_v(:,5)+xl_v(:,1));                                  yface5(:,4)= yl_v(:,1);                                                   zface5(:,4) = 0.5*(zl_v(:,5)+zl_v(:,1));
xcentrface5(:,1)= 0.25*( xface5(:,1)+ xface5(:,2)+ xface5(:,3)+ xface5(:,4)); ycentrface5(:,1)= yl_v(:,1);                                         zcentrface5(:,1) =  0.25*(zface5(:,1)+ zface5(:,2)+ zface5(:,3)+ zface5(:,4));

xface6(:,1)= 0.5*(xl_v(:,4)+xl_v(:,3));                                  yface6(:,1)= yl_v(:,7);                                                   zface6(:,1) = 0.5*(zl_v(:,4)+zl_v(:,3));
xface6(:,2)= 0.5*(xl_v(:,3)+xl_v(:,7));                                  yface6(:,2)= yl_v(:,7);                                                   zface6(:,2) = 0.5*(zl_v(:,3)+zl_v(:,7));
xface6(:,3)= 0.5*(xl_v(:,7)+xl_v(:,8));                                  yface6(:,3)= yl_v(:,7);                                                   zface6(:,3) = 0.5*(zl_v(:,7)+zl_v(:,8));
xface6(:,4)= 0.5*(xl_v(:,8)+xl_v(:,4));                                  yface6(:,4)= yl_v(:,7);                                                   zface6(:,4) = 0.5*(zl_v(:,8)+zl_v(:,4));
xcentrface6(:,1)= 0.25*( xface6(:,1)+ xface6(:,2)+ xface6(:,3)+ xface6(:,4)); ycentrface6(:,1)= yl_v(:,7);                                         zcentrface6(:,1) =  0.25*(zface6(:,1)+ zface6(:,2)+ zface6(:,3)+ zface6(:,4));

xmid(:) = 1./6*(xcentrface1(:,1)+xcentrface2(:,1)+xcentrface3(:,1)+xcentrface4(:,1)+xcentrface5(:,1)+xcentrface6(:,1));
ymid(:) = 1./6*(ycentrface1(:,1)+ycentrface2(:,1)+ycentrface3(:,1)+ycentrface4(:,1)+ycentrface5(:,1)+ycentrface6(:,1));
zmid(:) = 1./6*(zcentrface1(:,1)+zcentrface2(:,1)+zcentrface3(:,1)+zcentrface4(:,1)+zcentrface5(:,1)+zcentrface6(:,1));

% first
xl_s(:,1,1) = xl_v(:,1);          yl_s(:,1,1) = yl_v(:,1);         zl_s(:,1,1) = zl_v(:,1);
xl_s(:,1,2) = xface1(:,1);        yl_s(:,1,2) = yface1(:,1);       zl_s(:,1,2) = zface1(:,1);
xl_s(:,1,3) = xcentrface1(:,1);   yl_s(:,1,3) = ycentrface1(:,1);  zl_s(:,1,3) = zcentrface1(:,1);
xl_s(:,1,4) = xface1(:,4);        yl_s(:,1,4) = yface1(:,4);       zl_s(:,1,4) = zface1(:,4);
xl_s(:,1,5) = xface5(:,4);        yl_s(:,1,5) = yface5(:,4);       zl_s(:,1,5) = zface5(:,4);
xl_s(:,1,6) = xcentrface5(:,1);   yl_s(:,1,6) = ycentrface5(:,1);  zl_s(:,1,6) = zcentrface5(:,1);
xl_s(:,1,7) = xmid(:);            yl_s(:,1,7) = ymid(:);           zl_s(:,1,7) = zmid(:);
xl_s(:,1,8) = xcentrface4(:,1);   yl_s(:,1,8) = ycentrface4(:,1);  zl_s(:,1,8) = zcentrface4(:,1);
% second
xl_s(:,2,1) = xface1(:,1);        yl_s(:,2,1) = yface1(:,1);       zl_s(:,2,1) = zface1(:,1);
xl_s(:,2,2) = xl_v(:,2);          yl_s(:,2,2) = yl_v(:,2);         zl_s(:,2,2) = zl_v(:,2);
xl_s(:,2,3) = xface1(:,2);        yl_s(:,2,3) = yface1(:,2);       zl_s(:,2,3) = zface1(:,2);
xl_s(:,2,4) = xcentrface1(:,1);   yl_s(:,2,4) = ycentrface1(:,1);  zl_s(:,2,4) = zcentrface1(:,1);
xl_s(:,2,5) = xcentrface5(:,1);   yl_s(:,2,5) = ycentrface5(:,1);  zl_s(:,2,5) = zcentrface5(:,1);
xl_s(:,2,6) = xface2(:,4);        yl_s(:,2,6) = yface2(:,4);       zl_s(:,2,6) = zface2(:,4);
xl_s(:,2,7) = xcentrface2(:,1);   yl_s(:,2,7) = ycentrface2(:,1);  zl_s(:,2,7) = zcentrface2(:,1);
xl_s(:,2,8) = xmid(:);            yl_s(:,2,8) = ymid(:);           zl_s(:,2,8) = zmid(:);
% third
xl_s(:,3,1) = xcentrface1(:,1);   yl_s(:,3,1) = ycentrface1(:,1);  zl_s(:,3,1) = zcentrface1(:,1);
xl_s(:,3,2) = xface1(:,2);        yl_s(:,3,2) = yface1(:,2);       zl_s(:,3,2) = zface1(:,2);
xl_s(:,3,3) = xl_v(:,3);          yl_s(:,3,3) = yl_v(:,3);         zl_s(:,3,3) = zl_v(:,3);
xl_s(:,3,4) = xface1(:,3);        yl_s(:,3,4) = yface1(:,3);       zl_s(:,3,4) = zface1(:,3);
xl_s(:,3,5) = xmid(:);            yl_s(:,3,5) = ymid(:);           zl_s(:,3,5) = zmid(:);
xl_s(:,3,6) = xcentrface2(:,1);   yl_s(:,3,6) = ycentrface2(:,1);  zl_s(:,3,6) = zcentrface2(:,1);
xl_s(:,3,7) = xface6(:,2);        yl_s(:,3,7) = yface6(:,2);       zl_s(:,3,7) = zface6(:,2);
xl_s(:,3,8) = xcentrface6(:,1);   yl_s(:,3,8) = ycentrface6(:,1);  zl_s(:,3,8) = zcentrface6(:,1);
% fourth
xl_s(:,4,1) = xface4(:,1);        yl_s(:,4,1) = yface4(:,1);       zl_s(:,4,1) = zface4(:,1);
xl_s(:,4,2) = xcentrface1(:,1);   yl_s(:,4,2) = ycentrface1(:,1);  zl_s(:,4,2) = zcentrface1(:,1);
xl_s(:,4,3) = xface1(:,3);        yl_s(:,4,3) = yface1(:,3);       zl_s(:,4,3) = zface1(:,3);
xl_s(:,4,4) = xl_v(:,4);          yl_s(:,4,4) = yl_v(:,4);         zl_s(:,4,4) = zl_v(:,4);
xl_s(:,4,5) = xcentrface4(:,1);   yl_s(:,4,5) = ycentrface4(:,1);  zl_s(:,4,5) = zcentrface4(:,1);
xl_s(:,4,6) = xmid(:);            yl_s(:,4,6) = ymid(:);           zl_s(:,4,6) = zmid(:);
xl_s(:,4,7) = xcentrface6(:,1);   yl_s(:,4,7) = ycentrface6(:,1);  zl_s(:,4,7) = zcentrface6(:,1);
xl_s(:,4,8) = xface6(:,4);        yl_s(:,4,8) = yface6(:,4);       zl_s(:,4,8) = zface6(:,4);
% fifth
xl_s(:,5,1) = xface4(:,4);        yl_s(:,5,1) = yface4(:,4);       zl_s(:,5,1) = zface4(:,4);
xl_s(:,5,2) = xcentrface5(:,1);   yl_s(:,5,2) = ycentrface5(:,1);  zl_s(:,5,2) = zcentrface5(:,1);
xl_s(:,5,3) = xmid(:);            yl_s(:,5,3) = ymid(:);           zl_s(:,5,3) = zmid(:);
xl_s(:,5,4) = xcentrface4(:,1);   yl_s(:,5,4) = ycentrface4(:,1);  zl_s(:,5,4) = zcentrface4(:,1);
xl_s(:,5,5) = xl_v(:,5);          yl_s(:,5,5) = yl_v(:,5);         zl_s(:,5,5) = zl_v(:,5);
xl_s(:,5,6) = xface5(:,3);        yl_s(:,5,6) = yface5(:,3);       zl_s(:,5,6) = zface5(:,3);
xl_s(:,5,7) = xcentrface3(:,1);   yl_s(:,5,7) = ycentrface3(:,1);  zl_s(:,5,7) = zcentrface3(:,1);
xl_s(:,5,8) = xface4(:,3);        yl_s(:,5,8) = yface4(:,3);       zl_s(:,5,8) = zface4(:,3);
% sixth
xl_s(:,6,1) = xcentrface5(:,1);   yl_s(:,6,1) = ycentrface5(:,1);  zl_s(:,6,1) = zcentrface5(:,1);
xl_s(:,6,2) = xface2(:,4);        yl_s(:,6,2) = yface2(:,4);       zl_s(:,6,2) = zface2(:,4);
xl_s(:,6,3) = xcentrface2(:,1);   yl_s(:,6,3) = ycentrface2(:,1);  zl_s(:,6,3) = zcentrface2(:,1);
xl_s(:,6,4) = xmid(:);            yl_s(:,6,4) = ymid(:);           zl_s(:,6,4) = zmid(:);
xl_s(:,6,5) = xface5(:,3);        yl_s(:,6,5) = yface5(:,3);       zl_s(:,6,5) = zface5(:,3);
xl_s(:,6,6) = xl_v(:,6);          yl_s(:,6,6) = yl_v(:,6);         zl_s(:,6,6) = zl_v(:,6);
xl_s(:,6,7) = xface2(:,3);        yl_s(:,6,7) = yface2(:,3);       zl_s(:,6,7) = zface2(:,3);
xl_s(:,6,8) = xcentrface3(:,1);   yl_s(:,6,8) = ycentrface3(:,1);  zl_s(:,6,8) = zcentrface3(:,1);
% seventh
xl_s(:,7,1) = xmid(:);            yl_s(:,7,1) = ymid(:);           zl_s(:,7,1) = zmid(:);
xl_s(:,7,2) = xcentrface2(:,1);   yl_s(:,7,2) = ycentrface2(:,1);  zl_s(:,7,2) = zcentrface2(:,1);
xl_s(:,7,3) = xface2(:,2);        yl_s(:,7,3) = yface2(:,2);       zl_s(:,7,3) = zface2(:,2);
xl_s(:,7,4) = xcentrface6(:,1);   yl_s(:,7,4) = ycentrface6(:,1);  zl_s(:,7,4) = zcentrface6(:,1);
xl_s(:,7,5) = xcentrface3(:,1);   yl_s(:,7,5) = ycentrface3(:,1);  zl_s(:,7,5) = zcentrface3(:,1);
xl_s(:,7,6) = xface2(:,3);        yl_s(:,7,6) = yface2(:,3);       zl_s(:,7,6) = zface2(:,3);
xl_s(:,7,7) = xl_v(:,7);          yl_s(:,7,7) = yl_v(:,7);         zl_s(:,7,7) = zl_v(:,7);
xl_s(:,7,8) = xface6(:,3);        yl_s(:,7,8) = yface6(:,3);       zl_s(:,7,8) = zface6(:,3);
% eightth
xl_s(:,8,1) = xcentrface4(:,1);   yl_s(:,8,1) = ycentrface4(:,1);  zl_s(:,8,1) = zcentrface4(:,1);
xl_s(:,8,2) = xmid(:);            yl_s(:,8,2) = ymid(:);           zl_s(:,8,2) = zmid(:);
xl_s(:,8,3) = xcentrface6(:,1);   yl_s(:,8,3) = ycentrface6(:,1);  zl_s(:,8,3) = zcentrface6(:,1);
xl_s(:,8,4) = xface6(:,4);        yl_s(:,8,4) = yface6(:,4);       zl_s(:,8,4) = zface6(:,4);
xl_s(:,8,5) = xface4(:,3);        yl_s(:,8,5) = yface4(:,3);       zl_s(:,8,5) = zface4(:,3);
xl_s(:,8,6) = xcentrface3(:,1);   yl_s(:,8,6) = ycentrface3(:,1);  zl_s(:,8,6) = zcentrface3(:,1);
xl_s(:,8,7) = xface6(:,3);        yl_s(:,8,7) = yface6(:,3);       zl_s(:,8,7) = zface6(:,3);
xl_s(:,8,8) = xl_v(:,8);          yl_s(:,8,8) = yl_v(:,8);         zl_s(:,8,8) = zl_v(:,8);

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
        [jac,invjac_v,phi,dphidx_v,dphidy_v,dphidz_v] = deriv3D(sigpt,tigpt,ligpt,xl_m,yl_m,zl_m);
        rhs = gauss_source3D(sigpt,tigpt,ligpt,xl_m,yl_m,zl_m);
        
        fdem(:,:,subelt) = fdem(:,:,subelt) + wght*rhs(:).*phi.*jac(:);
        ael(:,:,:,subelt) =ael(:,:,:,subelt)+ wght*(dphidx_v.*permute(dphidx_v,[1 3 2]) + .....
            dphidy_v.*permute(dphidy_v,[1 3 2]) + .....
            dphidz_v.*permute(dphidz_v,[1 3 2])).*invjac_v(:);
        % end of Gauss point loop
    end
    
    % end of subdivided element loop
end
fdem = permute(fdem,[1,3,2]);
ael  = permute(ael,[1,4,2,3]);
% manual assembly of subelement contributions
% --------------------------------------------------- first face ------------------------

%cedroid, node number 6
ae(:,1,1)  = ael(:,1,3,3) + ael(:,2,4,4) + ael(:,3,1,1) + ael(:,4,2,2);
ae(:,1,2)  = ael(:,2,4,7) + ael(:,3,1,6);
ae(:,1,4)  = ael(:,4,2,5) + ael(:,1,3,8);
ae(:,1,5)  = ael(:,1,3,6) + ael(:,2,4,5);
ae(:,1,6)  = ael(:,3,1,8) + ael(:,4,2,7);
fde(:,1)   = fdem(:,1,3) + fdem(:,2,4) + fdem(:,3,1) + fdem(:,4,2);
% -------------------------------------------------- Second face ------------------------
% cedroid, node number 2
ae(:,2,1)  = ael(:,2,7,4) + ael(:,3,6,1);
ae(:,2,2)  = ael(:,2,7,7) + ael(:,6,3,3) + ael(:,3,6,6) + ael(:,7,2,2);
ae(:,2,3)  = ael(:,6,3,8) + ael(:,7,2,5);
ae(:,2,5)  = ael(:,2,7,5) + ael(:,6,3,1);
ae(:,2,6)  = ael(:,3,6,8) + ael(:,7,2,4);
ae(:,2,7)  = ael(:,2,7,8) + ael(:,3,6,5) + ael(:,6,3,4) + ael(:,7,2,1);
fde(:,2)   = fdem(:,2,7) + fdem(:,3,6) + fdem(:,7,2) + fdem(:,6,3);
% -------------------------------------------------- third face ------------------------

% centroid node number 10
% centroid node number 3
ae(:,3,2) = ael(:,6,8,3)+ael(:,7,5,2);
ae(:,3,3) = ael(:,5,7,7) + ael(:,6,8,8) + ael(:,7,5,5) + ael(:,8,6,6);
ae(:,3,4) = ael(:,5,7,4) + ael(:,8,6,1);
ae(:,3,5) = ael(:,5,7,2);
ae(:,3,6) = ael(:,8,6,3);
ae(:,3,7) = ael(:,5,7,3) + ael(:,6,8,4) + ael(:,7,5,1) + ael(:,8,6,2);
fde(:,3)  = fdem(:,5,7) + fdem(:,6,8) + fdem(:,7,5) + fdem(:,8,6);
% -------------------------------------------------- fourth face ------------------------
% centroid node, number 4
ae(:,4,1)  = ael(:,1,8,3) + ael(:,4,5,2);
ae(:,4,3) = ael(:,5,4,7) + ael(:,8,1,6);
ae(:,4,4) = ael(:,1,8,8) + ael(:,4,5,5) + ael(:,5,4,4) + ael(:,8,1,1);
ae(:,4,5)  = ael(:,1,8,6);
ae(:,4,6) = ael(:,4,5,7) + ael(:,8,1,3);
ae(:,4,7) = ael(:,1,8,7) + ael(:,4,5,6)+ael(:,5,4,3)+ael(:,8,1,2);
fde(:,4)    = fdem(:,1,8) + fdem(:,4,5) + fdem(:,8,1) + fdem(:,5,4);
% center nodes of fith sixth face and the centre node all 8 elements


% node nummber 5 (center node of the fifth face),   node 19 (center node
% of the sixth face) node 14 center node
% ----------------------------------------------------------------------

ae(:,5,1)  = ael(:,1,6,3) + ael(:,2,5,4);
ae(:,5,2)  = ael(:,2,5,7) + ael(:,6,1,3);
ae(:,5,3) = ael(:,5,2,7) + ael(:,6,1,8);
ae(:,5,4) = ael(:,1,6,8) + ael(:,5,2,4);
ae(:,5,5)  = ael(:,1,6,6) + ael(:,2,5,5) + ael(:,5,2,2) + ael(:,6,1,1);
ae(:,5,7) = ael(:,1,6,7) + ael(:,2,5,8) + ael(:,5,2,3) + ael(:,6,1,4);


ae(:,6,1)  = ael(:,4,7,2) + ael(:,3,8,1);
ae(:,6,2)  = ael(:,3,8,6) + ael(:,7,4,2);
ae(:,6,3) = ael(:,8,3,6)+ael(:,7,4,5);
ae(:,6,4) = ael(:,4,7,5) + ael(:,8,3,1);
ae(:,6,6) = ael(:,3,8,8) + ael(:,4,7,7) + ael(:,7,4,4) + ael(:,8,3,3);
ae(:,6,7) = ael(:,3,8,5) + ael(:,4,7,6) + ael(:,7,4,1) + ael(:,8,3,2);


ae(:,7,1)  = ael(:,1,7,3) + ael(:,2,8,4) + ael(:,3,5,1) + ael(:,4,6,2);
ae(:,7,2)  = ael(:,2,8,7) + ael(:,3,5,6) + ael(:,6,4,3) + ael(:,7,1,2);
ae(:,7,3) = ael(:,5,3,7) + ael(:,6,4,8) + ael(:,7,1,5) + ael(:,8,2,6);
ae(:,7,4) = ael(:,1,7,8) + ael(:,4,6,5) + ael(:,5,3,4) + ael(:,8,2,1);
ae(:,7,5)  = ael(:,1,7,6) + ael(:,2,8,5) + ael(:,5,3,2) + ael(:,6,4,1);
ae(:,7,6) = ael(:,3,5,8) + ael(:,4,6,7) + ael(:,7,1,4) + ael(8,2,3);
ae(:,7,7) = ael(:,1,7,7) + ael(:,2,8,8) + ael(:,3,5,5) + ael(:,4,6,6) + ael(:,5,3,3) + ael(:,6,4,4) + ael(:,7,1,1)+ael(:,8,2,2);

fde(:,5)   = fdem(:,1,6) + fdem(:,2,5) + fdem(:,6,1) + fdem(:,5,2);

fde(:,6)   = fdem(:,4,7) + fdem(:,3,8) + fdem(:,7,4) + fdem(:,8,3);

fde(:,7)   = fdem(:,1,7) + fdem(:,2,8) + fdem(:,3,5) + fdem(:,4,6) +....
    fdem(:,5,3) + fdem(:,6,4) + fdem(:,7,1)+ fdem(:,8,2);


% compute face residuals
res_face = faceres_Q1_Q1halfreduced3D(xyz,ev,ebound3D,q1sol3D,fcx,hx,hy,hz);
fprintf('Q1(h/2) estimate |  res_face = %7.4e\n', norm(res_face))
fe = fde - res_face;
% debug
%         if tout>0,
%         fprintf('right-hand side vector \n')
%         disp(squeeze(fe(elt,1:end)))
%         end

% solve for local estimate (sequential code)
%  for ielem = 1:nel
%      elerr(:,ielem) = squeeze(ae(ielem,1:nn,1:nn))\(fe(ielem,1:nn)');
%  end
%
% solve for local estimate (vectorized code)
% LDLT factorization
% redefine nn
nn=7;
dd=zeros(nel,nn); rr=zeros(nel,nn);
for kk=1:nn-1
    for pp=1:kk-1;
        rr(1:nel,pp)=dd(1:nel,pp).*ae(1:nel,kk,pp);
    end
    dd(1:nel,kk)=ae(1:nel,kk,kk);
    for pp=1:kk-1;
        dd(1:nel,kk)= dd(1:nel,kk)- ae(1:nel,kk,pp).*rr(1:nel,pp);
    end
    for ii=kk+1:nn
        for pp=1:kk-1;
            ae(1:nel,ii,kk)=ae(1:nel,ii,kk)-ae(1:nel,ii,pp).*rr(1:nel,pp);
        end
        ae(1:nel,ii,kk)=ae(1:nel,ii,kk)./dd(1:nel,kk);
    end
end
for pp=1:nn-1;
    rr(1:nel,pp)=dd(1:nel,pp).*ae(1:nel,nn,pp);
end
dd(1:nel,nn)=ae(1:nel,nn,nn);
for pp=1:nn-1;
    dd(1:nel,nn)= dd(1:nel,nn)- ae(1:nel,nn,pp).*rr(1:nel,pp);
end
% overwrite diagonal entries
for kk=1:nn
    ae(1:nel,kk,kk)= dd(1:nel,kk);
end
%  forward-backward substitutions ...
xx = element_lusolve(ae,fe);
elerr=xx';

%
for ivtx=1:7,
    err_sq_el(:) = err_sq_el(:) + fe(:,ivtx) .* elerr(ivtx,:)';
end
etime=toc;
fprintf('error estimation took %6.3e seconds\n',etime)

return