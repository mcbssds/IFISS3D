% CUBE_DIFF solves Poisson problem on unit cube domain
% IFISS scriptfile: DJS; 29 July 2022.
% Copyright (c) 2022 G. Papanikos, C.E. Powell, D.J. Silvester

%------ define geometry
pde=1; domain=1;
fprintf('\n\nGrid generation for cube domain\n')
nc=default('grid parameter: 3 for underlying 8x8x8 grid (default is 16x16x16)',4);
cube_domain(1,nc)
load cube_grid.mat

%------ set up matrices
qmethod=default('Q1/Q2 approximation 1/2? (default Q1)',1);
% reference grid switch
if grid_type==1  % uniform grid
    savesol=default('save results for reference 1/0 (yes/no)? (default no)',0);
else, savesol=0; end
if qmethod ==2,
    [x,y,z,xyz] = q2grid3D(x,y,z,xyz,mv,bound3D);
    [A,M,f] = femq2_diff3D(xyz,mv);
else
    ev = q1grid3D(xyz,mv,bound3D);
    [A,M,f] = femq1_diff3D(xyz,ev);
end

%------ apply boundary conditions
[Agal,fgal] = nonzerobc3D(A,f,xyz,bound3D);

%------ save resulting system
fprintf('system saved in cube_diff.mat ...\n')
gohome
cd datafiles
save -v7.3 cube_diff.mat qmethod Agal M fgal xyz x y z

lin_system_choice = default('Choose between a direct or an iterative solver 1/0 (direct/iterative) (default 1)',1);
tic
if lin_system_choice
    fprintf('solving linear system using direct solver... \n')
    x_it=Agal\fgal;
else
    fprintf('solving linear system using iterative solver...  ')
    solve_it
end
etoc=toc; fprintf('Galerkin system solved in %8.3e seconds\n\n',etoc)

%------ compute a posteriori error estimate and plot solution
if qmethod==1,
    diffpost3D
    errplot3D(x_it,error_tot,ev,xyz,x,y,z,99),
elseif qmethod==2,
    plot3Dsol(x_it,x,y,z,100)
end

if savesol == 1
    save cube_diff.mat x_it error_tot fcx hx hy hz -append
end

