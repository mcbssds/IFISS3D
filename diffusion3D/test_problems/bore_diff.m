% BORE_DIFF solves Poisson problem on a borehole domain
% IFISS scriptfile: DJS; 14 August 2022.
% Copyright (c) 2022 G. Papanikos, C.E. Powell, D.J. Silvester

% define geometry
pde=1; domain=3;
fprintf('\n\nGrid generation for borehole domain\n')
nc=default('Grid parameter? enter 2 for 48x48x32 grid (default 1)',1);
bore_domain(1,nc)
load bore_grid.mat
grid_type = 2;

% set up matrices
qmethod=1; %-------- Q2 option not tested
% reference grid switch
if grid_type==1  % uniform grid
    savesol=default('Save results for reference 1/0 (yes/no)? (default no)',0);
else, savesol=0; end
if qmethod ==2,
    [x,y,z,xyz] = q2grid3D(x,y,z,xyz,mv,bound3D);
    [A,M,f] = femq2_diff3D(xyz,mv);
else
    ev = q1grid3D(xyz,mv,bound3D);
    [A,M,f] = femq1_diff3D(xyz,ev);
end

% boundary conditions
[Agal,fgal] = nonzerobc3D(A,f,xyz,bound3D);

% save resulting system
fprintf('System saved in borehole_diff.mat ...\n')
gohome
cd datafiles
save borehole_diff.mat qmethod Agal M fgal xyz x y z

lin_system_choice = default('Choose between iterative or direct solver 0/1 (default 0)',0);
tic
if lin_system_choice
    fprintf('\nSolving linear system using direct solver... \n')
    x_it=Agal\fgal;
else
    fprintf('\nSolving linear system using iterative solver...  ')
    solve_it
end
etoc=toc; fprintf('Galerkin system solved in %8.3e seconds\n\n',etoc)

% compute a posteriori error estimate (Q1 only) and plot solution
if qmethod==1,
    diffpost3D
    errplot3D_borehole(x_it,error_tot,ev,xyz,x,y,z,bd,99),
elseif qmethod==2,
    plot3Dsol_borehole(x_it,xyz,x,y,z,bd,100)
end

if savesol == 1
    save borehole_diff.mat x_it error_tot fcx hx hy hz -append
end

