% CUBE_DIFF solves Poisson problem on cube domain
% IFISS scriptfile: DJS; 25 September 2022.
% Copyright (c) 2022 G. Papanikos, C.E. Powell, D.J. Silvester

% define geometry
pde=1; domain=1;
fprintf('\n\nGrid generation for cube domain\n')
nc=default('Grid parameter: 3 for underlying 8x8x8 grid (default is 16x16x16)',4);
gtic=tic; cube_domain(1,nc), toc(gtic)
load cube_grid.mat

% set up matrices
qmethod=default('Q1/Q2 approximation 1/2? (default Q1)',1);
% reference grid switch
if grid_type==1  % uniform grid
    savesol=default('Save results for reference 1/0 (yes/no)? (default no)',0);
else, savesol=0; end
atic=tic;
if qmethod ==2,
    [x,y,z,xyz] = q2grid3D(x,y,z,xyz);
    [A,M,f] = femq2_diff3D(xyz,mv);
else
    ev = q1grid3D(xyz,mv);
    [A,M,f] = femq1_diff3D(xyz,ev);
end

% apply boundary conditions
[Agal,fgal] = nonzerobc3D(A,f,xyz,bound3D);
atoc=toc(gtic);
fprintf('Galerkin system assembled in %8.3e seconds\n\n',atoc)

lin_system_choice = default('Choose between direct or iterative solver 1/0 (direct/iterative) (default 1)',1);
stic=tic;
if lin_system_choice
    fprintf('\nSolving linear system using direct solver... \n')
    x_it=Agal\fgal;
else
    fprintf('\nSolving linear system using iterative solver...  ')
    solve_it
end
etoc=toc(stic); fprintf('Galerkin system solved in %8.3e seconds\n\n',etoc)

% compute a posteriori error estimate (Q1 only) and plot solution
if qmethod==1,
etic=tic; diffpost3D, rtoc=toc(etic);
    fprintf('Error estimated in %8.3e seconds\n\n',rtoc)
    errplot3D(x_it,error_tot,ev,xyz,x,y,z,99),
elseif qmethod==2,
    plot3Dsol(x_it,x,y,z,100)
end

if savesol == 1
fprintf('System saved in cube_diff.mat ...\n')
gohome
cd datafiles
save cube_diff.mat qmethod Agal M fgal xyz x y z x_it
end
