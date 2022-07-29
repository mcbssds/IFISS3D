%DIFFPOST3D driver for a posteriori error postprocessing
% This implements the full/reduced Q2 full/reduced Q1(h/2) error approximation scheme of Qifeng Liao
% IFISS scriptfile: DJS; 18 July 2022.
% Copyright (c) 2022 G.Papanikos, C.E. Powell, D.J. Silvester

fprintf('choose a specific error estimation space');
fprintf('\n     1.  Q2(h)')
fprintf('\n     2.  reduced Q2(h)')
fprintf('\n     3.  Q1(h/2)')
fprintf('\n     4.  reduced Q1(h/2)\n')

approx_space = default('default is reduced Q2 approximation',2);

fprintf('FAST a posteriori error estimation for Q1 element \n')
[hx,hy,hz,fcx,ebound3D] = facegen(xyz,ev,domain);
if  approx_space==1
    [errorsq_ele,elerr,fe,ae] = diffpostq1_q2_3D(xyz,ev,ebound3D,x_it,fcx,hx,hy,hz);
    error_tot = sqrt(errorsq_ele); errorest=norm(error_tot(:),2);
    fprintf('simplistic estimate energy error is %10.4e \n',errorest)
    
    if sn~=5
        % include the boundary correction
        [errorsq_cbc] = diffpost_bc3D(errorsq_ele,fe,xyz,ev,ebound3D);
        error_tot = sqrt(errorsq_cbc); errorest=norm(error_tot,2);
        fprintf('corrected estimate energy error is %10.4e \n',errorest)
    end
elseif approx_space==2
    [errorsq_ele,elerr,fe,ae] = diffpostq1_3D_q2_reduced(xyz,ev,ebound3D,x_it,fcx,hx,hy,hz);
    error_tot = sqrt(errorsq_ele); errorest=norm(error_tot(:),2);
    fprintf('simplistic estimate energy error is %10.4e \n',errorest)
    if sn~=5
        % include the boundary correction
        [errorsq_cbc] = diffpost_bc3D_q2_reduced(errorsq_ele,fe,xyz,ev,ebound3D);
        error_tot = sqrt(errorsq_cbc); errorest=norm(error_tot,2);
        fprintf('corrected estimate energy error is %10.4e \n',errorest)
    end
elseif approx_space==3
    [errorsq_ele,elerr,fe,ae,xl_m,yl_m,zl_m,xl_s,yl_s,zl_s] = diffpost_Q1_Q1half3D(xyz,ev,ebound3D,x_it,fcx,hx,hy,hz);
    error_tot = sqrt(errorsq_ele); errorest=norm(error_tot(:),2);
    fprintf('simplistic estimate energy error is %10.4e \n',errorest)
    if sn~=5
        % include the boundary correction
        [errorsq_cbc] = diffpost_bc3D_Q1half(errorsq_ele,fe,xyz,ev,ebound3D,xl_m,yl_m,zl_m,xl_s,yl_s,zl_s);
        error_tot = sqrt(errorsq_cbc); errorest=norm(error_tot,2);
        fprintf('corrected estimate energy error is %10.4e \n',errorest)
    end
elseif approx_space == 4
    [errorsq_ele,elerr,fe,ae,xl_m,yl_m,zl_m,xl_s,yl_s,zl_s] = diffpost_Q1_Q1halfreduced3D(xyz,ev,ebound3D,x_it,fcx,hx,hy,hz);
    error_tot = sqrt(errorsq_ele); errorest=norm(error_tot(:),2);
    fprintf('simplistic estimate energy error is %10.4e \n',errorest)
    if sn~=5
        % include the boundary correction
        [errorsq_cbc] = diffpost_bc3D_Q1halfreduced(errorsq_ele,fe,xyz,ev,ebound3D,xl_m,yl_m,zl_m,xl_s,yl_s,zl_s);
        error_tot = sqrt(errorsq_cbc); errorest=norm(error_tot,2);
        fprintf('corrected estimate energy error is %10.4e \n',errorest)
    end
else
    fprintf('Invalid input!');
end

